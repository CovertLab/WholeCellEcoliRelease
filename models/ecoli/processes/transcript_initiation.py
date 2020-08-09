# -*- coding: utf-8 -*-
"""
    @author: John Mason
    @organization: Covert Lab, Department of Bioengineering, Stanford University
    @date: Created 4/26/14
    LastEditors: Hwrn
    LastEditTime: 2020-08-09 11:44:11
    ImportPath: models.ecoli.processes.transcript_initiation
    Description:
        TranscriptInitiation
        Transcription initiation sub-model.
    TODO:
        - use transcription units instead of single genes
        - match sigma factors to promoters
        - implement transcriptional regulation
        - modulate initiation probabilities as a function of gene copy number
        - match measured levels of active RNA polymerase instead of initiating to completion
"""

import numpy as np
import scipy.sparse

import wholecell.processes.process
from wholecell.utils import units

import itertools

class TranscriptInitiation(wholecell.processes.process.Process):
    """ TranscriptInitiation """

    _name = "TranscriptInitiation"

    # Constructor
    def __init__(self):
        super(TranscriptInitiation, self).__init__()

    def initialize(self, sim, sim_data):
        super(TranscriptInitiation, self).initialize(sim, sim_data)

        # Load parameters
        self.fracActiveRnapDict = sim_data.process.transcription.rnapFractionActiveDict
        self.rnaLengths = sim_data.process.transcription.rnaData["length"]
        self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict

        # 转录招募矩阵 Create recruitment matrix for transcription regulation
        recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
        recruitmentData = sim_data.process.transcription_regulation.recruitmentData
        # 根据 "hI" 和 "hJ" 值将 "hV" 值放入矩阵
        self.recruitmentMatrix = scipy.sparse.csr_matrix(
                (recruitmentData["hV"], (recruitmentData["hI"], recruitmentData["hJ"])),
                shape = recruitmentData["shape"]
            )

        # 基因扰动, 用于计算合成概率 Determine changes from genetic perturbations
        self.genetic_perturbations = {}
        perturbations = getattr(sim_data, "genetic_perturbations", {})
        if len(perturbations) > 0:
            rnaIdxs, synthProbs = zip(*[(int(np.where(sim_data.process.transcription.rnaData["id"] == rnaId)[0]), synthProb) for rnaId, synthProb in sim_data.genetic_perturbations.iteritems()])
            fixedSynthProbs = [synthProb for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
            fixedRnaIdxs = [rnaIdx for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
            self.genetic_perturbations = {"fixedRnaIdxs": fixedRnaIdxs, "fixedSynthProbs": fixedSynthProbs}

        # If initiationShuffleIdxs does not exist, set value to None
        self.shuffleIdxs = getattr(sim_data.process.transcription, 'initiationShuffleIdxs', None)

        # Views
        self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
        self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
        self.chromosomes = self.bulkMoleculeView('CHROM_FULL[c]')
        self.recruitmentView = self.bulkMoleculesView(recruitmentColNames)

        # 都是数组, 长度相同 | ID Groups
        self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
        self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
        self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']
        self.isRRna = sim_data.process.transcription.rnaData['isRRna']
        self.isMRna = sim_data.process.transcription.rnaData["isMRna"]
        self.isTRna = sim_data.process.transcription.rnaData["isTRna"]
        self.isRProtein = sim_data.process.transcription.rnaData['isRProtein']
        self.isRnap = sim_data.process.transcription.rnaData['isRnap']
        self.isRegulated = np.array([
            1 if x[:-3] in sim_data.process.transcription_regulation.targetTf or x in perturbations else 0
            for x in sim_data.process.transcription.rnaData["id"]
        ], dtype = np.bool)
        # 或操作, 记录重要 RNA
        self.setIdxs = self.isRRna | self.isTRna | self.isRProtein | self.isRnap | self.isRegulated

        # 基因合成概率 Synthesis probabilities for different categories of genes
        self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction  # 所以多加一个 s 是什么意思 @Hwrn: ?
        self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
        self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

    def calculateRequest(self):
        # Get all inactive RNA polymerases
        self.inactiveRnaPolys.requestAll()

        # 根据转录调节计算合成概率 Calculate synthesis probabilities based on transcription regulation
        # np.array rnaSynthProb: $\vec{p_{synth}}$ 向量
        self.rnaSynthProb = self.recruitmentMatrix.dot(self.recruitmentView.total())
        if self.genetic_perturbations:  # not {}
            self.rnaSynthProb[self.genetic_perturbations["fixedRnaIdxs"]] = self.genetic_perturbations["fixedSynthProbs"]
        regProbs = self.rnaSynthProb[self.isRegulated]

        # 最低值为 0.0 Adjust probabilities to not be negative
        self.rnaSynthProb[self.rnaSynthProb < 0] = 0.0
        self.rnaSynthProb /= self.rnaSynthProb.sum()
        # Warning: useless assert
        assert np.all(self.rnaSynthProb >= 0), "Have negative RNA synthesis probabilities"

        # 根据环境选择参数 Adjust synthesis probabilities depending on environment
        current_nutrients = self._external_states['Environment'].nutrients

        # synthProbFractions: $f_{act}$. 与环境有关
        synthProbFractions = self.rnaSynthProbFractions[current_nutrients]
        # 使用 $\tag{22}$ 计算合成概率
        self.rnaSynthProb[self.isMRna] *= synthProbFractions["mRna"] / self.rnaSynthProb[self.isMRna].sum()
        self.rnaSynthProb[self.isTRna] *= synthProbFractions["tRna"] / self.rnaSynthProb[self.isTRna].sum()
        self.rnaSynthProb[self.isRRna] *= synthProbFractions["rRna"] / self.rnaSynthProb[self.isRRna].sum()
        self.rnaSynthProb[self.isRegulated] = regProbs
        self.rnaSynthProb[self.isRProtein] = self.rnaSynthProbRProtein[current_nutrients]
        self.rnaSynthProb[self.isRnap] = self.rnaSynthProbRnaPolymerase[current_nutrients]
        # 概率均 >= 0
        self.rnaSynthProb[self.rnaSynthProb < 0] = 0
        # 对于相对不重要的 RNA, 合成概率是剩余的加权平分
        scaleTheRestBy = (1. - self.rnaSynthProb[self.setIdxs].sum()) / self.rnaSynthProb[~self.setIdxs].sum()
        self.rnaSynthProb[~self.setIdxs] *= scaleTheRestBy

        # 改变初始速率 Shuffle initiation rates if we're running the variant that calls this
        # (In general, this should lead to a cell which does not grow and divide)
        if self.shuffleIdxs is not None:
            self.rnaSynthProb = self.rnaSynthProb[self.shuffleIdxs]

        # int fracActiveRnap: $f_{act}$ INPUT
        self.fracActiveRnap = self.fracActiveRnapDict[current_nutrients]
        # 延伸率
        self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[current_nutrients]

    def evolveState(self):
        self.writeToListener("RnaSynthProb", "rnaSynthProb", self.rnaSynthProb)

        # 无基因组则返回 no synthesis if no chromosome
        if self.chromosomes.total()[0] == 0:
            return

        # 根据概率计算激活的转录酶数 Calculate RNA polymerases to activate based on probabilities
        # $p_{act}$
        self.activationProb = self._calculateActivationProb(
            self.fracActiveRnap, self.rnaLengths, self.rnaPolymeraseElongationRate, self.rnaSynthProb
        )
        # rnaPolyToActivate: $c_{RNAP, b} = p_{act} \cdot c_{RNA, f}$
        rnaPolyToActivate = np.int64(self.activationProb * self.inactiveRnaPolys.count())
        if rnaPolyToActivate == 0:
            return

        #### 生长控制 Growth control code ####
        # 对合成概率的多项式分布进行抽样, 允许一个基因多次转录 Sample a multinomial distribution of synthesis probabilities to determine what molecules are initialized
        # TODO 时间粒度可以进一步细化, 直到同一个位点只能被一个蛋白占据. 是否有必要这么仔细?
        # self.randomState: np.random.RandomState
        # \vec{n_{init}} = \text{multinomial}(c_{RNAP,b}, \vec{p_{synth}})
        nNewRnas = self.randomState.multinomial(rnaPolyToActivate, self.rnaSynthProb)

        # Build list of RNA indexes
        # 每位对应一个开始转录的基因
        rnaIndexes = np.empty(rnaPolyToActivate, np.int64)
        startIndex = 0
        # 仅对非零值进行操作
        nonzeroCount = (nNewRnas > 0)
        # np.arange(nNewRnas.size)[nonzeroCount] 提取非零值的指针
        for rnaIndex, counts in itertools.zip(np.arange(nNewRnas.size)[nonzeroCount], nNewRnas[nonzeroCount]):
            rnaIndexes[startIndex:startIndex+counts] = rnaIndex
            startIndex += counts

        # 激活酶 Create the active RNA polymerases
        activeRnaPolys = self.activeRnaPolys.moleculesNew("activeRnaPoly", rnaPolyToActivate)
        activeRnaPolys.attrIs(rnaIndex=rnaIndexes)
        # 更新可用的转录酶的数量
        self.inactiveRnaPolys.countDec(nNewRnas.sum())

        # Write outputs to listeners
        self.writeToListener("RibosomeData", "rrn16S_produced", nNewRnas[self.is_16SrRNA].sum())
        self.writeToListener("RibosomeData", "rrn23S_produced", nNewRnas[self.is_23SrRNA].sum())
        self.writeToListener("RibosomeData", "rrn5S_produced", nNewRnas[self.is_5SrRNA].sum())

        self.writeToListener("RibosomeData", "rrn16S_init_prob", nNewRnas[self.is_16SrRNA].sum() / float(nNewRnas.sum()))
        self.writeToListener("RibosomeData", "rrn23S_init_prob", nNewRnas[self.is_23SrRNA].sum() / float(nNewRnas.sum()))
        self.writeToListener("RibosomeData", "rrn5S_init_prob", nNewRnas[self.is_5SrRNA].sum() / float(nNewRnas.sum()))

        self.writeToListener("RibosomeData", "total_rna_init", nNewRnas.sum())

        self.writeToListener("RnapData", "didInitialize", nNewRnas.sum())
        self.writeToListener("RnapData", "rnaInitEvent", nNewRnas)

    def _calculateActivationProb(self, fracActiveRnap, rnaLengths, rnaPolymeraseElongationRate, synthProb):
        """
            description:
                计算激活率
                p_{act} = \frac{r \cdot f_{act}}{1 - f_{act}}
            @param fracActiveRnap $f_{act}$, 激活转录酶的比例
            @param rnaLengths $\vec{}$
            @param rnaPolymeraseElongationRate $\vec{}$
                以上两者用于计算转录所需时间
            @param synthProb RNA 合成概率
            @return <int> 激活率, 根据 1.1.6.1 计算
                范围: [0, 1]
                Edge case: 100% RNA polymerase activity (relevant to paper investigations).
        """
        if fracActiveRnap == 1:
            return 1.

        # **计算终止期待值 expectedTerminationRate** | Calculate expected RNAP termination rate based on RNAP elongation rate
        # \vec{转录所需时间} = \vec{长度}/转录延伸率 | Vector of times required to transcribe each transcript
        allTranscriptionTimes = 1. / rnaPolymeraseElongationRate * rnaLengths
        # 时间片 = self.timeStepSec() * units.s | Vector of numbers of timesteps required to transcribe each transcript
        timesteps = (1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber()
        # \vec{转录所需时间片数量} 向上舍入整数
        allTranscriptionTimestepCounts = np.ceil(timesteps)
        # # 概率加权后的平均转录所需时间片 (期待值) | Average number of timesteps required to transcribe a transcript, weighted by synthesis probabilities
        # averageTranscriptionTimestepCounts = np.dot(synthProb, allTranscriptionTimestepCounts)
        # 终止期待值 $r$ | Average number of terminations in one timestep for one transcript
        expectedTerminationRate = 1. / np.dot(synthProb, allTranscriptionTimestepCounts)

        # **对转录中的酶, 考虑早期终止** (因此会有更多的可用的转录酶) | Modify given fraction of active RNAPs to take into account early terminations in between timesteps
        # 此概率随 (待转录 ?) 长度降低而增加, 在倒数第二个周期时仍 < 0.5 | Vector of probabilities an "active" RNAP will in effect be "inactive" because it has terminated during a timestep
        allFractionTimeInactive = 1 - timesteps / allTranscriptionTimestepCounts
        # 同样进行加权 | Average probability of an "active" RNAP being in effect "inactive", weighted by synthesis probabilities
        averageFractionTimeInactive = np.dot(synthProb, allFractionTimeInactive)
        # 更高的有效起始率 | 可能 >1 ? | New higher "goal" for fraction of active RNAP, considering that the "effective" fraction is lower than what the listener sees
        effectiveFracActiveRnap = fracActiveRnap * 1 / (1 - averageFractionTimeInactive)

        # $p_{act} = \frac{r \cdot f_{act}}{1 - f_{act}}$ | Return activation probability that will balance out the expected termination rate
        return effectiveFracActiveRnap * expectedTerminationRate / (1 - effectiveFracActiveRnap)
