#!/usr/bin/env -S bash -x
TEST_NAME=test_python_obs

#before running this script, set the following environment variables:
# TEST_ROOT_DIR [root directory of the project]
[ -z "${TEST_ROOT_DIR}" ] && { echo "environment variable TEST_ROOT_DIR not defined"; exit 1; }
# TEST_WORK_DIR [working directory to generate temporary output files]
[ -z "${TEST_WORK_DIR}" ] && { echo "environment variable TEST_WORK_DIR not defined"; exit 1; }
# TEST_SCRIPT_DIR [directory where test scripts are stored]
[ -z "${TEST_SCRIPT_DIR}" ] && { echo "environment variable TEST_SCRIPT_DIR not defined"; exit 1; }
# TEST_EXECUTABLE [path to the executable to generate output]
[ -z "${TEST_EXECUTABLE}" ] && { echo "environment variable TEST_EXECUTABLE not defined"; exit 1; }

#init variables
TEST_PYTHONPATH="${TEST_ROOT_DIR}/opt/python"
TEST_EVAL="python ${TEST_SCRIPT_DIR}/assets/test_python_obs.py ${TEST_PYTHONPATH}"

#write task file
cat > ${TEST_WORK_DIR}/${TEST_NAME}.xml <<- EOM
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
        <lattice name="honeycomb" range="3"/>
        <model name="honeycomb-heisenberg" symmetry="SU2">
            <j>1.0</j>
        </model>
    </parameters>
    <measurements>
        <measurement name="correlation"/>
    </measurements>
</task>
EOM

function cleanup {
    for EXT in xml obs ldf checkpoint data ; do
        rm -f ${TEST_WORK_DIR}/${TEST_NAME}.${EXT}
    done
}

#run executable
${TEST_EXECUTABLE} -f ${TEST_WORK_DIR}/${TEST_NAME}.xml

#evaluate test
trap 'cleanup ; exit 1' ERR
${TEST_EVAL} ${TEST_WORK_DIR}/${TEST_NAME}.obs

#cleanup
cleanup