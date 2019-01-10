filename=$1

workDir="/project/etp4/mherrmann/imageEvaluation/"

root -l -q -x ${workDir}'referenceNparse.C("'${workDir}'ROpanels/'${filename}'.txt","'${workDir}'zeroImages.txt","'${workDir}'parsed")'

root -l -q -x ${workDir}'combineCMMnPi.C("'${workDir}'parsed/'${filename}'_dif.txt","/project/etpdaq/cmm/camLeftPinsNmarker.txt")'

root -q -x ${workDir}'globalFit.C("'${workDir}'cmmPoints/'${filename}'_dif_cmm.txt")'

cd plotter

root -q -x ${workDir}'plotter/alignmentPlotter.C("'${workDir}'residuals/'${filename}'_dif_cmm_residuals.txt")'
