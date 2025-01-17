#!/bin/bash

# unpack PGAP reference data
tar xzf /srv/input-2022-04-14.build6021.tgz

ls -lh

# gzip -dv "${5}"

# cwltool -h
cwltool --timestamps --disable-color --preserve-entire-environment --outdir /srv/output /pgap/pgap/pgap.cwl "/srv/${2}"

# rm "${6}"

ANNOT=$(find . -name annot.gff)

echo ${ANNOT}
echo ${4}

mv "${ANNOT}" "${4}"

# mv /svr/output/annot.gff "${4}"

ls -lh
