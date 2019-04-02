#!/bin/sh
#
# rep-unit-typing3.3.sh
#

#
_TTL_PROGRAM1=RepUnitTyping3.3.py
_REPTYPE_DIR=$HOME/RepUnitTyping
_REP_OUT=RepUnit_out
#
        echo
        echo "\033[1;33m<<< $_TTL_PROGRAM1 >>>\033[m"
        echo
        echo "*** When .fastq or .fastq.gz files are targets, $_TTL_PROGRAM1 runs on file names containing "_R1", "_1_" or "_1.fastq" first, and then "_R2", "_2_" or "_2.fastq" second"
        echo "*** When .fasta (or .fa) files are targets, this rule is not applied. Single fasta files are read one after another"
        echo
        echo "Source file type: [ fq / fastq / fq.gz / fastq.gz / fa / fna / fasta ] ?"
	echo "(ex. type: fastq.gz )"
	echo -n "type: "
	read _FASTQ_QGZ_EXT
        echo
        echo Source $_FASTQ_QGZ_EXT file directory:
        echo "(ex. $HOME/MTB_fastq_gz )"
        read _FASTQ_QGZFULLPATH
        echo
        echo ALL $_FASTQ_QGZ_EXT files from: $_FASTQ_QGZFULLPATH
        echo
	cd $_FASTQ_QGZFULLPATH
        ls ./*.$_FASTQ_QGZ_EXT
        echo
        echo -n "No. files: "
        ls ./*.$_FASTQ_QGZ_EXT | wc -l
	echo
        echo
	date +%T
        _CURR_DATE=`date +%Y%m%d`
        #
        cd $_REPTYPE_DIR

        echo "*** ref subdirectory: "
        ls ./ref
        echo
        echo "Choose query filename (ex. rep_unit.fasta ) in the ref subdirectory: "
        read _QUNIT
        echo
        echo "Output directory name under the $_REPTYPE_DIR directory: "
        echo "(ex. RepUnit_out )"
        read _REP_OUT
        if [ -e $_REP_OUT ] ; then
            echo "Output files are saved in the $_REP_OUT subdirectory"
        else
            echo "Output files are saved in a new $_REP_OUT subdirectory"
            mkdir "$_REP_OUT"
        fi
        echo
        echo -n "Prediction mode? [y/n] "
        while :
        do
            read _ANS1
            case "$_ANS1" in
            y | Y )
            _PRED="-p"
            echo "Prediction mode!"
            break
            ;;
            n | N )
            _PRED=""
            echo "Non-prediction mode!"
            break
            ;;
            * )
            echo "Please enter y or n!"
            ;;
            esac
        done

        echo
        echo -n "Usually NOT necessary but detail mode? [y/n] "
        while :
        do
            read _ANS1
            case "$_ANS1" in
            y | Y )
            _DETAIL="-d"
            echo "Detail mode!"
            break
            ;;
            n | N )
            _DETAIL=""
            echo "Non-detail mode!"
            break
            ;;
            * )
            echo "Please enter y or n!"
            ;;
            esac
        done



        while :
        do
	#cd 
cat <<EOT

        *Input file extension: $_FASTQ_QGZ_EXT
        *Query fasta file: $_QUNIT
        *Option(s): $_PRED $_DETAIL 
        
        ### Select the running mode

        [1] swift mode (for fq, fastq, fq.gz, or fastq.gz input files)

        [2] non-swift mode (for fq, fastq, fq.gz, or fastq.gz files)
       
        [3] fasta mode (for fa, fna, or fasta files)

        [q] quit without analysis


EOT
        echo -n "Input Number: " 

                read _ANS_MENU
                case "$_ANS_MENU" in
                1)
		_NUM=0
        	#cd $_REPTYPE_DIR
        	echo "$_TTL_PROGRAM1 runs under the $_REPTYPE_DIR directory."
                echo
                for _FILE1 in $_FASTQ_QGZFULLPATH/*$_FASTQ_QGZ_EXT
                do
                        if [ `echo "$_FILE1" | grep "_R1"` ] || [ `echo "$_FILE1" | grep "_1_"` ] || [ `echo "$_FILE1" | grep "_1.fastq"` ]
                        then
                                _NUM=`expr $_NUM + 1`
                                echo +++++++++++++++++
                                echo "[$_NUM]"
				echo "  <<<<<< $_FILE1 >>>>>>"
				_FILE2=`echo "$_FILE1" | sed -e "s%_R1%_R2%" | sed -e "s%_1_%_2_%" | sed -e "s%_1.fastq%_2.fastq%"`
				echo "  <<<<<< $_FILE2 >>>>>>"
				date +%T
                                echo +++++++++++++++++
                                echo
#####################################
                               # output files contain _R1, _1_ or_1.fastq
				echo "Typing results are in the $_REP_OUT directory"
				_FOUT="${_CURR_DATE}RepUnitTyping"
				_COMMAND="python $_TTL_PROGRAM1 $_FILE1 $_FILE2 -q $_QUNIT -O $_REP_OUT -o $_FOUT $_PRED $_DETAIL"
				echo "command: $_COMMAND"
				$_COMMAND

				echo "<< swift mode: $_FILE1 and $_FILE2 => ${_REP_OUT}/$_FOUT file >> done! "
#~/RepUnitTyping$ python RepUnitTyping.py $HOME/MTB_fastq/DN-049_S6_L001_R1_001.trimmed.fastq $HOME/MTB_fastq/DN-049_S6_L001_R2_001.trimmed.fastq -o DN-049_S6.RepUnitTyping

#####################################
                        fi
                done

                break
                ;;

                2)
		echo "<< non-swift mode >>"
                _NUM=0
                #cd $_REPTYPE_DIR
                echo "$_TTL_PROGRAM1 runs under the $_REPTYPE_DIR directory."
                echo
                for _FILE1 in $_FASTQ_QGZFULLPATH/*$_FASTQ_QGZ_EXT
                do
                        if [ `echo "$_FILE1" | grep "_R1"` ] || [ `echo "$_FILE1" | grep "_1_"` ] || [ `echo "$_FILE1" | grep "_1.fastq"` ]
                        then
                                _NUM=`expr $_NUM + 1`
                                echo +++++++++++++++++
                                echo "[$_NUM]"
                                echo "  <<<<<< $_FILE1 >>>>>>"
                                _FILE2=`echo "$_FILE1" | sed -e "s%_R1%_R2%" | sed -e "s%_1_%_2_%" | sed -e "s%_1.fastq%_2.fastq%"`
                                echo "  <<<<<< $_FILE2 >>>>>>"
                                date +%T
                                echo +++++++++++++++++
                                echo
#####################################
                                #_R1 or _1_ is contained
                                echo "Typing results are in the $_REP_OUT directory"
                                _FOUT="${_CURR_DATE}RepUnitTyping"
                                _COMMAND="python $_TTL_PROGRAM1 -s off $_FILE1 $_FILE2 -q $_QUNIT -O $_REP_OUT -o $_FOUT $_PRED $_DETAIL"
                                echo "command: $_COMMAND"
                                $_COMMAND

                                echo "<< non-swift mode: $_FILE1 and $_FILE2 => ${_REP_OUT}/$_FOUT file >> done! "
#~/RepUnitTyping$ python RepUnitTyping.py -s off /media/nkrit/HDPC-UT/0928-Malawi-again/ERR037476_1.fastq.gz /media/nkrit/HDPC-UT/0928-Malawi-again/ERR037476_1.fastq.gz -o 20160929RepUnitTyping

#####################################
                        fi
                done

                break
                ;;

                3)
                _NUM=0
                #cd $_REPTYPE_DIR
                echo "$_TTL_PROGRAM1 runs under the $_REPTYPE_DIR directory."
                echo "fasta file extension should be .$_FASTQ_QGZ_EXT !!"
                echo
                for _FILE1 in $_FASTQ_QGZFULLPATH/*$_FASTQ_QGZ_EXT
                do
                                _NUM=`expr $_NUM + 1`
                                echo +++++++++++++++++
                                echo "[$_NUM]"
				echo "  <<<<<< $_FILE1 >>>>>>"
				date +%T
                                echo +++++++++++++++++
                                echo
#####################################
				echo "Typing results are in the $_REP_OUT directory"
				_FOUT="${_CURR_DATE}RepUnitTyping"
				_COMMAND="python $_TTL_PROGRAM1 $_FILE1 -q $_QUNIT -O $_REP_OUT -o $_FOUT --seq $_PRED $_DETAIL"
				echo "command: $_COMMAND"
				$_COMMAND

				echo "<< fasta mode: $_FILE1 => ${_REP_OUT}/$_FOUT file >> done! "
#~/RepUnitTyping$ python RepUnitTyping.py $HOME/MTB_fasta/DN-049.fasta -o DN-049_S6.RepUnitTyping --seq -d

#####################################
                done

                break
                ;;

                q | Q)
                #cd $_REPTYPE_DIR
                echo [End]
                break
                ;;

                *)
                echo Please enter 1, 2 or q!
                ;;
                esac
        done

		echo
                echo All files in ${_REPTYPE_DIR}/${_REP_OUT}:
                ls ./${_REP_OUT}
		cd
