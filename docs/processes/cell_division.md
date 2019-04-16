\paragraph{Cell division}

\subparagraph{Model implementation.}
Cell division is modeled in the \texttt{ChromosomeElongation} process and \texttt{CellDivision} listener in the \emph{E. coli} model. A Helmstetter-Cooper type model of chromosome replication initiation is coupled to cell division, inspired by work from Wallden \emph{et al.}  \cite{Wallden2016}. Chromosome replication initiation occurs at a fixed mass per origin of replication. Each initiation event is coupled to a cell division event after a constant period of time consisting of one round of chromosome replication and cytokinesis. Importantly, this constant period of time can span multiple cell division events.

Cell division itself is modeled as a binomial process where each daughter cell has an equal probability of inheriting the contents of the mother cell. The exception to this is if two chromosomes are present before cell division---each daughter is guaranteed to get one.

Algorithms \ref{dna_replication_termination_algorithm} and \ref{cell_division_algorithm} provide implementation details.\\

\begin{algorithm}[H]
\caption{Algorithm for cell division}
\label{cell_division_algorithm}
\SetKwInOut{Input}{Input}\SetKwInOut{Result}{Result}
\SetKwFunction{peek}{peek}
\SetKwFunction{pop}{pop}
\SetKwFunction{rand}{rand}
\SetKwFunction{mod}{mod}
\SetKwFunction{floor}{floor}
\SetKwFunction{randint}{randint}

  \Input{$d_{queue}$ a double ended queue data structure that stores time(s) cell division should be triggered}
     \Input{$c_i$ counts of all molecules in simulation at cell division where $i = 1$ \KwTo $n_{species}$}
     \Input{$p$ binomial partition coefficient}
     \Input{$n_{chrom}$ number of chromosome molecules}
     \Input{\rand{} returns a random number from a uniform distribution between 0 and 1}
     \Input{\randint{} returns a random integer either 0 or 1}

  \If{$t > d_{queue}$.\peek{}}{
      \textbf{1.} Trigger division and remove division time.\\
        \-\hspace{1cm} $d_{queue}$.\pop{}\\
      \textbf{2.} Divide bulk contents of cell binomially. Number partitioned into daughter one is stored in $n_{daughter,1}$ and to daughter two in $n_{daughter,2}.$\\
          \For{$i = 1$ \KwTo $n_{species}$}{
              $n_{daughter,1} = 0$\\
              \For{$j = 1$ \KwTo $c_i$}{
                  \If{\rand{} $> p$}{
                      $n_{daughter,1} = n_{daughter,1} + 1$
                  }
              }
              $n_{daughter,2} = c_i - n_{daughter,1}$
          }
        \textbf{3.} Divide chromosome in binary manner. All replication forks and origins of replication associated with a chromosome molecule are partitioned as well. Number of chromosome molecules partitioned into daughter one is stored in $n_{chrom,daughter,1}$ and to daughter two in $n_{chrom,daughter,2}.$\\
          \eIf{\mod{$n_{chrom}$,2}}{
              $n_{chrom,daughter,1} = \frac{n_{chrom}}{2}$\\
            }
            {
              $n_{chrom,daughter,1} = $\floor{$\frac{n_{chrom}}{2}$} + \randint{}\\
            }
            $n_{chrom,daughter,2} = n_{chrom} - n_{chrom,daughter,1}$

    }

\Result{Cell division is triggered at C+D time after DNA replication initiation. Contents of mother cell is divided between two daughter cells conserving mass.}
\end{algorithm}

\newpage
%\subparagraph*{Associated files} 
\textbf{Associated files}

\begin{table}[h!]
 \centering
 \begin{tabular}{c c c} 
 \hline
 \texttt{wcEcoli} Path & File & Type \\
 \hline
\texttt{wcEcoli/models/ecoli/processes} & \texttt{replication\_elongation.py} & process \\
\texttt{wcEcoli/models/ecoli/listeners} & \texttt{cell\_division.py} & listener \\
 \hline
\end{tabular}
\caption[Table of files for transcription regulation]{Table of files for transcription regulation.}
\end{table}


\subparagraph{Difference from \emph{M. genitalium} model.}
The \emph{E. coli} model is not yet a gene complete model and many of the mechanistic details of cell division are not implemented as they were in the \emph{M. genitalium} model. In \emph{E. coli}, cytokinesis, septation, and chromosome segregation are all not modeled explicitly.  However, cell division in our \textit{E. coli} model is consistent with growth at multiple growth rates, which was not the case in \textit{M. genitalium}.
