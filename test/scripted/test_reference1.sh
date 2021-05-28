#!/usr/bin/env bash
TEST_NAME=test_reference1

#before running this script, set the following environment variables:
# TEST_WORK_DIR [working directory to generate temporary output files]
[ -z "${TEST_WORK_DIR}" ] && { echo "environment variable TEST_WORK_DIR not defined"; exit 1; }
# TEST_SCRIPT_DIR [directory where test scripts are stored]
[ -z "${TEST_SCRIPT_DIR}" ] && { echo "environment variable TEST_SCRIPT_DIR not defined"; exit 1; }
# TEST_EXECUTABLE [path to the executable to generate output]
[ -z "${TEST_EXECUTABLE}" ] && { echo "environment variable TEST_EXECUTABLE not defined"; exit 1; }

#init variables
TEST_EVAL="python ${TEST_SCRIPT_DIR}/assets/test_eval.py"
TEST_REF_DIR="${TEST_SCRIPT_DIR}/assets"

#write task files
for CORE in SU2 XYZ TRI ; do 
cat > ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.xml <<- EOM
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
        <lattice name="hyperkagome" range="3"/>
        <model name="hyperkagome-heisenberg" symmetry="${CORE}">
            <j>1.0</j>
        </model>
    </parameters>
    <measurements>
        <measurement name="correlation" />
    </measurements>
</task>
EOM
done

function cleanup {
    for CORE in SU2 XYZ TRI ; do
        for EXT in xml obs ldf checkpoint data ; do
            rm -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${EXT}
        done
    done
}

#run executable
for CORE in SU2 XYZ TRI ; do 
    ${TEST_EXECUTABLE} -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.xml
done

#evaluate test
trap 'cleanup ; exit 1' ERR
${TEST_EVAL} FILE ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs ${TEST_REF_DIR}/${TEST_NAME}.ref
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorDD ${TEST_WORK_DIR}/${TEST_NAME}.XYZ.obs XYZCorDD
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.XYZ.obs XYZCorXX
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.XYZ.obs XYZCorYY
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.XYZ.obs XYZCorZZ
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorDD ${TEST_WORK_DIR}/${TEST_NAME}.TRI.obs TRICorDD
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.TRI.obs TRICorXX
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.TRI.obs TRICorYY
${TEST_EVAL} OBJECT ${TEST_WORK_DIR}/${TEST_NAME}.SU2.obs SU2CorZZ ${TEST_WORK_DIR}/${TEST_NAME}.TRI.obs TRICorZZ

#cleanup
cleanup
