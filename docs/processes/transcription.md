
\paragraph{Transcription}

\subparagraph{Model implementation.}
Transcription occurs through the action of two processes in the model: \texttt{TranscriptInitiation} and \texttt{TrancriptElongation}. \texttt{TranscriptInitiation} models the binding of RNA polymerase  to each gene. The number of initiation events per gene is proportional to the number of free RNA polymerases weighted by each gene's synthesis probability. Details are in Algorithm \ref{transcript_initiation_algorithm}.

\texttt{TranscriptElongation} models nucleotide polymerization into RNA molecules by RNA polymerases. Polymerization occurs across all polymerases simultaneously and resources are allocated to maximize the progress of all polymerases up to the limit of the expected polymerase elongation rate and available nucleotides. The termination of RNA elongation occurs once a RNA polymerase has reached the end of the annotated gene. Details are in Algorithm \ref{transcript_elongation_algorithm}.


\subparagraph{Difference from \emph{M. genitalium} model.}
The \emph{M. genitalium} model modeled RNA polymerase as existing in 4 states: free, non-specifically bound on a chromosome, bound to a promoter, and actively transcribing a gene. The \emph{E. coli} model simplifies this by assuming RNA polymerase exists in two states: free and actively transcribing. Every time step, free RNA polymerase transitions to the actively transcribing state to maintain an experimentally-observed active fraction of RNA polymerase. The \emph{E. coli} model does not yet include sigma, elongation or termination factors. The \textit{E. coli} model also currently treats each gene as its own transcription unit.\\

\begin{algorithm}[H]
\caption{Algorithm for RNA polymerase initiation on DNA}
\label{transcript_initiation_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}
\SetKwFunction{min}{min}
\SetKwFunction{multinomial}{multinomial}
  \Input{$f_{act}$ fraction of RNA polymerases that are active}
    \Input{$r$ expected termination rate for active RNA polymerases}
  \Input{$v_{\text{synth,i}}$ RNA synthesis probability for each gene where $i = 1$ \KwTo $n_{gene}$}
    \Input{$c_{RNAP,f}$ count of free RNA polymerase}
    \Input{\multinomial{} function that draws samples from a multinomial distribution}
    
  \textbf{1.} Calculate probability ($p_{act}$) of a free RNA polymerase binding to a gene.\\
    \-\hspace{1cm} $p_{act} = \frac{f_{act} \cdot r}{1 - f_{act}}$
    
    \textbf{2.} Calculate the number of RNA polymerases that will bind and activate ($c_{RNAP,b}$).\\
    \-\hspace{1cm} $c_{RNAP,b} = p_{act} \cdot c_{RNAP,f}$
    
    \textbf{3} Sample multinomial distribution $c_{RNAP,b}$ times weighted by $v_{synth,i}$ to determine which genes receive a RNA polymerase and initiate ($n_{init,i}$).\\
    \-\hspace{1cm} $n_{init,i} =$ \multinomial{$c_{RNAP,b}, v_{\text{synth}, i}$}
    
    \textbf{4} Assign $n_{init,i}$ RNA polymerases to gene $i$. Decrement free RNA polymerase counts.\\

    \Result{RNA polymerases bind to genes based on the number of free RNA polymerases and the synthesis probability for each gene.}

\end{algorithm}
\newpage



\begin{algorithm}[H]
\caption{Algorithm for mRNA elongation and termination}
\label{transcript_elongation_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}
\SetKwFunction{min}{min}
\SetKwFunction{all}{all}

  \Input{$e$ expected RNA polymerase elongation rate in given environment}
  \Input{$L_i$ length of each gene $i = 1$ \KwTo $n_{gene}$ for each coding gene.}
    \Input{$p_j$ gene position of RNA polymerase $j = 1$ \KwTo $n_{RNAP}$}
    \Input{$c_{nuc,k}$ counts of nucleotide $k = 1$ \KwTo $4$}
  \Input{$\delta t$ length of current time step}
    
    \tcc{Elongate RNA transcripts up to limits of sequence or nucleotides}
    \For{each RNA polymerase j on gene i}{
      \textbf{1.} Based on RNA polymerase position $p_j$ on a gene $i$ and maximal elongation rate $e$ determine ``stop condition'' ($s_j$) for RNA polymerase $j$ assuming no nucleotide limitation.\\
      \-\hspace{1cm} $s_j = $ \min{$p_j + e \cdot \delta t$, $L_i$}\\
        Stop condition is either maximal elongation rate scaled by the time step or the full length of sequence (i.e. the RNA polymerase will terminate in this time step).\\
    
      \textbf{2.} Derive sequence between RNA polymerase position ($p_j$) and stop condition ($s_j$).

  \textbf{3.} Based on derived sequence calculate the number of nucleotides required to polymerize sequence $c^{req}_{nuc,k}$.\\
    
    \textbf{4.} Elongate up to limits:\\
    \eIf{\all{$c^{req}_{nuc,k} < c_{nuc,k}$}}{
      Update the position of each polymerase to stop position\\
      \-\hspace{1cm} $p_j = s_j$
    }
    {
      \textbf{4a.} Attempt to elongate all RNA fragments.\\
        \textbf{4b.} Update position of each polymerase to maximal position given the limitation of $c_{nuc,k}$.
    }

    \textbf{5.} Update counts of $c_{nuc,k}$ to reflect polymerization usage.
    }
    \tcc{Terminate RNA polymerases that have reached the end of their gene}
    \For{each RNA polymerase j on gene i}{
      \If{$p_j$ == $L_i$}{
          \textbf{1.} Increment count of RNA that corresponds to elongating RNA transcript that has terminated.\\
            
            \textbf{2.} Increment free RNA polymerase counts.
        
        }
    }
    
    \Result{Each RNA transcript is elongated up to the limit of available gene sequence, expected elongation rate, or nucleotide limitation. RNA polymerases that reach the end of their genes are terminated and released.}
\end{algorithm}

\newpage
%\subparagraph*{Associated data}
\textbf{Associated data}

\begin{table}[h!]
 \centering
 \label{transcript_initiation_table}
 \begin{tabular}{c c c c c} 
 \hline
 Parameter & Symbol & Units & Value & Reference \\
 \hline
  Active fraction of RNAP & $f_{act}$ & - & 0.20 (growth-dependent) & \cite{bremer2008modulation} \\
  RNA synthesis probability$^{(1)}$ & $p_{synth}$ & - & [0, 0.015] & See Table \ref{files_transcription} \\
  RNAP elongation rate & $e$ & nt/s & 50 (growth-dependent) & \cite{bremer2008modulation} \\
 \hline
\end{tabular}
\caption[Table of parameters for Transcript Initiation and Elongation]{Table of parameters for Transcript Initiation and Elongation processes.\\$^{(1)}$RNA synthesis probabilities were calculated as the relative fraction of RNA production (which is equal to the RNA degradation) for a given gene.}
\end{table}

%\subparagraph*{Associated files}
\textbf{Associated files}

\begin{table}[h!]
 \centering
 \scriptsize
 \begin{tabular}{c c c} 
 \hline
 \texttt{wcEcoli} Path & File & Type \\
 \hline
\texttt{wcEcoli/models/ecoli/processes} & \texttt{transcript\_initiation.py} & process \\
\texttt{wcEcoli/models/ecoli/processes} & \texttt{transcript\_elongation.py} & process \\
\texttt{wcEcoli/reconstruction/ecoli/dataclasses/process} & \texttt{transcription.py} & data \\
 \hline
\end{tabular}
\caption[Table of files for transcription]{Table of files for transcription.}
\label{files_transcription}
\end{table}
