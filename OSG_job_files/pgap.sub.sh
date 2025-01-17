universe = vanilla
executable = PGAPWrapper.sh
log = PGAPLogFiles/out.$(v01).log
output = PGAPLogFiles/out.$(v01).out
error = PGAPLogFiles/out.a.err

# OSG has different requirement syntax than CHTC
# Kernel version comes from major.minor.patch in the format:
# major * 10000 + minor * 1000 + patch
# heron comes from r-base which uses debian:latest, which
# as of 20190722 is "buster" from 4.19.105
# as per advice, 31000 will keep jobs away from RHEL 6

requirements = Arch == "X86_64" && HAS_SINGULARITY == True && OSG_HOST_KERNEL_VERSION >= 31000 && SINGULARITY_CAN_USE_SIF == False
request_cpus = 6
request_memory = 64GB
request_disk = 64GB

# may need to change this
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/ncbi/pgap:2022-04-14.build6021"

# IF using a docker container, only send the RScript and any associated data files
transfer_input_files = YAMLFiles/$(v02), \
                        YAMLFiles/$(v03), \
                        $(v05), \
			stash:///osgconnect/public/npcooley/NCBI/input-2022-04-14.build6021.tgz

max_materialize = 4000

arguments = $(v01) $(v02) $(v03) $(v04) $(v05) $(v06)

queue v01, v02, v03, v04, v05, v06 from PGAPJobMap.txt

