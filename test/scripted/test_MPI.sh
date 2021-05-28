#!/usr/bin/env bash
TEST_NAME=test_MPI

#before running this script, set the following environment variables:
# TEST_WORK_DIR [working directory to generate temporary output files]
[ -z "${TEST_WORK_DIR}" ] && { echo "environment variable TEST_WORK_DIR not defined"; exit 1; }
# TEST_SCRIPT_DIR [directory where test scripts are stored]
[ -z "${TEST_SCRIPT_DIR}" ] && { echo "environment variable TEST_SCRIPT_DIR not defined"; exit 1; }
# TEST_EXECUTABLE [path to the executable to generate output]
[ -z "${TEST_EXECUTABLE}" ] && { echo "environment variable TEST_EXECUTABLE not defined"; exit 1; }

# TEST_MPIEXEC_EXECUTABLE [path to the mpi wrapper]
[ -z "${TEST_MPIEXEC_EXECUTABLE}" ] && { echo "environment variable TEST_MPIEXEC_EXECUTABLE not defined"; exit 1; }
# TEST_MPIEXEC_NUMPROC_FLAG [path to the mpi wrapper flag to specify number of ranks]
[ -z "${TEST_MPIEXEC_NUMPROC_FLAG}" ] && { echo "environment variable TEST_MPIEXEC_NUMPROC_FLAG not defined"; exit 1; }

#init variables
TEST_EVAL="python ${TEST_SCRIPT_DIR}/assets/test_eval.py"

#write task files
for CORE in SU2 XYZ TRI ; do 
    for MODE in MPI NMPI ; do 
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
</task>
EOM
    done
done

function cleanup {
    for CORE in SU2 XYZ TRI ; do
        for MODE in MPI NMPI ; do 
            for EXT in xml obs ldf checkpoint data ; do
                rm -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.${MODE}.${EXT}
            done
        done
    done
}

#run executable
for CORE in SU2 XYZ TRI ; do 
    ${TEST_EXECUTABLE} -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.NMPI.xml
    ${TEST_MPIEXEC_EXECUTABLE} ${TEST_MPIEXEC_NUMPROC_FLAG} 4 --use-hwthread-cpus ${TEST_EXECUTABLE} -f ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.MPI.xml
done

#evaluate test
trap 'cleanup ; exit 1' ERR
for CORE in SU2 XYZ TRI ; do 
    ${TEST_EVAL} FILE ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.NMPI.obs ${TEST_WORK_DIR}/${TEST_NAME}.${CORE}.MPI.obs
done

#cleanup
cleanup