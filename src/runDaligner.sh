#!/bin/bash
dir=/ebio/abt6_projects7/small_projects/mdubarry
fileReference=allSequencesNodes.fasta
fileRead=tmpRead28.fasta

cd $dir/Documents/SampleProgram/bin/output
rm *.dam .d* dbReference.dbRead.las #rm all files
# Create Database
	#Reference
export PATH=$dir/software/DAZZ_DB:$PATH
fasta2DAM ./dbReference $dir/Documents/pbsim/$fileReference
	#Read
fasta2DAM ./dbRead $dir/Documents/pbsim/$fileRead

# Align read against reference
export PATH=$dir/software/DALIGNER:$PATH
HPCmapper $dir/Documents/SampleProgram/bin/output/dbReference $dir/Documents/SampleProgram/bin/output/dbRead -v | csh -v
#LAshow to look at the alignment
LAshow dbReference.dam dbRead.dam dbReference.dbRead.las -ca
# Convert .las to .sam using proovread
export PATH=$dir/software/proovread/bin:$PATH
export PERL5LIB=$dir/software/proovread/lib:$PERL5LIB
dazz2sam --las dbReference.dbRead.las --ref dbReference --qry dbRead > referenceRead.sam
#rm dbReference.dbRead.las #no need to keep this file, convert in sam format
