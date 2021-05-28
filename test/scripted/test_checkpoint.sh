#!/usr/bin/env bash
TEST_NAME=test_checkpoint

#before running this script, set the following environment variables:
# TEST_WORK_DIR [working directory to generate temporary output files]
[ -z "${TEST_WORK_DIR}" ] && { echo "environment variable TEST_WORK_DIR not defined"; exit 1; }
# TEST_SCRIPT_DIR [directory where test scripts are stored]
[ -z "${TEST_SCRIPT_DIR}" ] && { echo "environment variable TEST_SCRIPT_DIR not defined"; exit 1; }
# TEST_EXECUTABLE [path to the executable to generate output]
[ -z "${TEST_EXECUTABLE}" ] && { echo "environment variable TEST_EXECUTABLE not defined"; exit 1; }

#init variables
TEST_EVAL="python ${TEST_SCRIPT_DIR}/assets/test_eval.py"

#write task files
for CORE in SU2 XYZ TRI ; do 
    for MODE in CHKPNT NOCHKPNT ; do 
        if [ ${MODE} == CHKPNT ] ; then
            CUTOFF_MIN=0.5
        else
            CUTOFF_MIN=0.3
        fi
        cat > ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.xml <<- EOM
<?xml version="1.0" encoding="utf-8"?>
<task>
    <parameters>
        <frequency discretization="manual">
            <value>0.31812</value>
            <value>0.36329</value>
            <value>0.41812</value>
            <value>0.46329</value>
            <value>0.51334</value>
            <value>0.56880</value>
            <value>0.63024</value>
            <value>0.69833</value>
            <value>0.77378</value>
            <value>0.85737</value>
            <value>0.95</value>
            <value>1.0</value>
            <value>3.0</value>
            <value>10.0</value>
        </frequency>
        <cutoff discretization="exponential">
            <max>10</max>
            <min>${CUTOFF_MIN}</min>
            <step>0.9</step>
        </cutoff>
        <lattice name="triangular" range="3"/>
        <model name="triangular-heisenberg" symmetry="${CORE}">
            <j>1.0</j>
        </model>
    </parameters>
    <measurements>
        <measurement name="correlation" />
    </measurements>
</task>
EOM
    done
done

function cleanup {
    for CORE in SU2 XYZ TRI ; do
        for MODE in CHKPNT NOCHKPNT ; do 
            for EXT in xml obs ldf checkpoint data ; do
                rm -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.${EXT}
            done
        done
    done
}

#run executable
for CORE in SU2 XYZ TRI ; do 
    for MODE in CHKPNT NOCHKPNT ; do 
        ${TEST_EXECUTABLE} -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.xml
    done
done

#rewrite task file for resuming checkpoint
for CORE in SU2 XYZ TRI ; do 
    for MODE in CHKPNT ; do 
        cat > ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.xml <<- EOM
<?xml version="1.0" encoding="utf-8"?>
<task>
    <parameters>
        <frequency discretization="manual">
            <value>0.31812</value>
            <value>0.36329</value>
            <value>0.41812</value>
            <value>0.46329</value>
            <value>0.51334</value>
            <value>0.56880</value>
            <value>0.63024</value>
            <value>0.69833</value>
            <value>0.77378</value>
            <value>0.85737</value>
            <value>0.95</value>
            <value>1.0</value>
            <value>3.0</value>
            <value>10.0</value>
        </frequency>
        <cutoff discretization="exponential">
            <max>10</max>
            <min>0.3</min>
            <step>0.9</step>
        </cutoff>
        <lattice name="triangular" range="3"/>
        <model name="triangular-heisenberg" symmetry="${CORE}">
            <j>1.0</j>
        </model>
    </parameters>
    <measurements>
        <measurement name="correlation" />
    </measurements>
    <calculation status="running" startTime="1970-Jan-01 00:00:00" checkpointTime="1970-Jan-01 00:00:00" />
</task>
EOM
    done
done

#run executable
for CORE in SU2 XYZ TRI ; do 
    for MODE in CHKPNT ; do 
        ${TEST_EXECUTABLE} ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.xml
    done
done

#evaluate test
trap 'cleanup ; exit 1' ERR
for CORE in SU2 XYZ TRI ; do 
    ${TEST_EVAL} FILE ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.CHKPNT.obs ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.NOCHKPNT.obs
    ${TEST_EVAL} FILE ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.CHKPNT.obs ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.NOCHKPNT.obs
    ${TEST_EVAL} FILE ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.CHKPNT.obs ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.NOCHKPNT.obs
done

#cleanup
cleanup