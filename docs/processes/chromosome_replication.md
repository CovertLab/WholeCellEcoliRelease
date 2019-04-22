\paragraph{Chromosome replication}

\subparagraph{Model implementation.}
Chromosome replication occurs through three steps that are implemented in the \texttt{ChromosomeFormation} and \texttt{ChromosomeElongation} processes. First, a round of replication is initiated at a fixed cell mass per origin of replication and generally occurs once per cell cycle (see Algorithm~\ref{replication_init_algorithm}). Second, replication forks are elongated up to the maximal expected elongation rate, dNTP resource limitations, and template strand sequence  (see Algorithm~\ref{dna_replication_elongation_algorithm}). Finally, replication forks terminate once they reach the end of their template strand and the chromosome immediately decatinates forming two separate chromosome molecules  (see Algorithm~\ref{dna_replication_termination_algorithm}).\\

\begin{algorithm}[H]
\caption{Algorithm for DNA replication initiation}
\label{replication_init_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}

  \Input{$m_{cell}$ cell mass}
  \Input{$m_{critical}$ critical initiation mass}
  \Input{$n_{origin}$ number of origins of replication}
  \Input{$n_{fork,f}$ number of replication forks on forward strand}
  \Input{$n_{fork,r}$ number of replication forks on reverse strand}
  \Input{$n_{chromosome}$ number of chromosome molecules}
  \Input{$C$ length of C period}
  \Input{$D$ length of D period}
  
  \If{$\frac{m_{cell}}{n_{origin}} > m_{critical}$}{
    \eIf{$n_{origin} > 1$}{
      $n_{origin} = n_{origin} + \frac{n_{fork,f} + n_{fork,r}}{2} \cdot n_{chromosome}$
    }
    {
      $n_{origin} = n_{origin} + n_{chromosome}$
    }
    
    $n_{fork,f} = n_{fork,f} + n_{fork,f} \cdot n_{chromosome}$\\
    $n_{fork,r} = n_{fork,r} + n_{fork,r} \cdot n_{chromosome}$
  }
  
  \Result{When cell mass is larger than critical initiation mass $m_c$ another round of replication is initiated with correct number of replication forks}
\end{algorithm}

\newpage

\begin{algorithm}[H]
\caption{Algorithm for DNA replication elongation}
\label{dna_replication_elongation_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}
\SetKwFunction{min}{min}
\SetKwFunction{all}{all}

  \Input{$e$ maximal elongation rate of replication fork}
    \Input{$p_i$ position of forks on chromosome where $i = 1$ \KwTo $n_{fork}$}
  \Input{$\delta t$ length of current time step}
    \Input{$c_{dNTP,j}$ counts of dNTP where $j = 1$ \KwTo $4$ for dCTP, dGTP, dATP, dTTP}
  \Input{$L_k$ total length of each strand of chromosome from origin to terminus where $k = 1$ \KwTo $4$ for forward/complement and reverse/complement.}
    \For{each replication fork $i$ on sequence $k$}{
      \textbf{1.} Based on replication fork position $p_i$ and maximal elongation rate $e$ determine ``stop condition'' ($s_i$) for replication fork assuming no dNTP limitation.\\
      \-\hspace{1cm} $s_i = $ \min{$p_i + e \cdot \delta t$, $L_k$}\\
        Stop condition is either maximal elongation rate scaled by the time step or the full length of sequence (i.e. the fork will terminate in this time step).\\
    
      \textbf{2.} Derive sequence between replication fork position ($p_i$) and stop condition ($s_i$).

  \textbf{3.} Based on derived sequence calculate the number of dNTPs required to polymerize sequence $c^{req}_{dNTP,i}$.\\
    
    \textbf{4.} Elongate up to limits:\\
    \eIf{\all{$c^{req}_{dNTP,i} < c_{dNTP,j}$}}{
      Update the position of each replication fork to stop position\\
      \-\hspace{1cm} $p_i = s_i$
    }
    {
      Attempt to equally elongate each replication fork update position of each fork to maximal position given the limitation of $c_{dNTP,j}$.
    }

    \textbf{5.} Update counts of $c_{dNTP,j}$ to reflect polymerization usage.
    }
    
    \Result{Each replication fork is elongated up to the limit of available sequence, elongation rate, or dNTP limitation}
\end{algorithm}
\vspace{1cm}
\begin{algorithm}[H]
\caption{Algorithm for DNA replication termination}
\label{dna_replication_termination_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}
\SetKwFunction{push}{push}

  \Input{$p_i$ position of forks on chromosome where $i = 1$ \KwTo $n_{fork}$}
  \Input{$L_k$ total length of each strand of chromosome from origin to terminus where $k = 1$ \KwTo $4$ for forward/complement and reverse/complement}
    \Input{$d_{queue}$ a double ended queue data structure that stores time(s) cell division should be triggered}
    \Input{$D$ D-period of cell cycle (time between completion of chromosome replication and cell division)}
    \Input{$t$ Current simulation time}

  \For{each replication fork i on strand k}{
      \If{$p_i == L_k$}{
      \textbf{1.} Delete replication fork\\
      \textbf{2.} Divide remaining replication forks and origins of replication appropriately across the two new chromosome molecules\\
      \textbf{3.} Calculate time cell should trigger division based on current time of chromosome termination and push onto queue data structure\\
            \-\hspace{1cm} $d_{queue}$.\push{$t + D$}
      }
    }
\Result{Replication forks that have terminated are removed. A new chromosome molecule is created separating all remaining replication forks. Timer for D-period is started.}
\end{algorithm}

\vspace{1 cm}
%\subparagraph*{Associated files}
\textbf{Associated files}

\begin{table}[h!]
 \centering
 \scriptsize
 \begin{tabular}{c c c} 
 \hline
 \texttt{wcEcoli} Path & File & Type \\
 \hline
\texttt{wcEcoli/models/ecoli/processes} & \texttt{chromosome\_formation.py} & process \\
\texttt{wcEcoli/models/ecoli/processes} & \texttt{chromosome\_elongation.py} & process \\
\texttt{wcEcoli/reconstruction/ecoli/dataclasses/process} & \texttt{replication.py} & data \\
 \hline
\end{tabular}
\caption[Table of files for chromosome replication]{Table of files for chromosome replication.}
\end{table}


\subparagraph{Difference from \emph{M. genitalium} model.}
The physiology modeled is significantly different from what was implemented in the \emph{M. genitalium} model. Initiation of DNA replication in \emph{E. coli} no longer uses a DnaA based mechanistic model but instead uses a phenomenological model based on a constant mass per origin of replication triggering DNA replication initiation. The action of topoisomerases are not explicitly modeled. Replication forks no longer take into account all of the enzymes in the replisome but are point objects that traverse the chromosome sequence. Some differences exist because the \emph{E. coli} model is not yet a gene complete model.  More importantly, certain changes enabled significant modeling advances in the \emph{E. coli} model. These include modeling the DNA replication cycle over multiple growth rates, cell sizes, and conditions using a single unified framework, and enabling multiple rounds of replication to proceed simultaneously over multiple generations. Both advances were critical to the findings in this publication.


%\subparagraph*{Associated data}
\textbf{Associated data}

\begin{table}[h!]
 \centering
 \begin{tabular}{p{1.75in} p{1in} p{1in} p{1in} p{1in}} 
 \hline
 Parameter & Symbol & Units & Value & Reference \\
 \hline
Chromosome sequence & - & - & - & \cite{Blattner:1997wl} \\
Replication fork elongation rate & $e$ & nt/s & 967 & \cite{Bremer:1996uj} \\
Mass per origin at DNA replication initiation & $m_{critical}$ & origin/fg & [600,975] & Semi-quantitative fit \cite{Donachie:1968vp} \\
C period & $C$ & min & 40 & \cite{Neidhardt:1990tna} \\
D period & $D$ & min & 20 & \cite{Neidhardt:1990tna} \\

 \hline
\end{tabular}
\caption[Table of parameters for chromosome replication]{Table of parameters for chromosome replication process.
}
\end{table}
