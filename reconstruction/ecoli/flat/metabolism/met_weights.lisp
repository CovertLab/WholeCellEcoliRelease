; TODO: units on MW
; TODO: generalize the paths?
; TODO: instructions on how to call from command-line
(with-open-file (stream "~/sherlock/scratch/users/jmason42/wcEcoli/reconstruction/ecoli/flat/metabolites.tsv" ; open the file and give it a handle
		:direction :output ; specify that this is an output file, default is to read
		:if-exists :supersede ; overwrite any existing file
		)
	(format stream "\"~A\"~C\"~A\"~%" ; write the header (array)(character)(array)(newline)
		"Compound ID"
		#\tab ; no format string for tab, have to write the character
		"MW"
		)
	(loop for x in (get-class-all-instances '|Compounds|) ; loop over compounds
		do (format stream "\"~A\"~C~A~%" ; write out to the stream
			(get-frame-name x) ; get the name of the compound
			#\tab
			(get-slot-value x 'MOLECULAR-WEIGHT) ; get the molecular weight
			)
		)
	)

; example call
; $ ./pathway-tools -lisp
; : (load "~/sherlock/scratch/users/jmason42/wcEcoli/reconstruction/ecoli/flat/metabolism/met_weights.lisp")
