rule align_mtb:
	input:
		u_1 = OUTPUT + '{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
		u_2 = OUTPUT + '{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq'
	output:
		bam = OUTPUT +  '{project}/sample_output{sample}/{sample}_mtb.bam',
		btw_report = OUTPUT + '{project}/sample_output/{sample}/{sample}_mtb.bt2.txt'
	threads: 12
	params:
		MTB_REF = '/workdir/ac2763/cfRNA_pipeline/takara_human_V2/references/MTb.H37Rv.fasta'
	shell:
		"""
		bwa-mem -x {params.MTB_REF} \
			-1 {input.u_1} -2 {input.u_2} \
			-p {threads} -q \
			--dovetail 2> {output.btw_report} | samtools view -f3 -bh - > {output.bam};
		"""


rule align_mtb_star:
	input:
		u_1 = OUTPUT + '{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
                u_2 = OUTPUT + '{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq'
	output: bam = OUTPUT + '{project}/sample_output/{sample}/{sample}_MTB_Aligned.sortedByCoord.out.bam'
	params:
		prefix = OUTPUT + "{project}/sample_output/{sample}/{sample}_MTB_",
		STAR = MPATH + config['STAR'],
		STAR_IDX = '/workdir/ac2763/cfRNA_pipeline/takara_human_V2/references/star'
	threads: 10
	shell:
		"""
		{params.STAR} \
			--genomeDir {params.STAR_IDX} \
			--readFilesIn {input.u_1} {input.u_2} \
			--outFileNamePrefix {params.prefix} \
			--runThreadN {threads} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outSAMmultNmax 20 \
			--outSAMprimaryFlag AllBestScore \
			--outFilterMismatchNmax 999 \
			--outFilterMismatchNoverReadLmax 0.04 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000
		"""
rule count_pos:
	input: OUTPUT + '{project}/sample_output/{sample}/{sample}_MTB_Aligned.sortedByCoord.out.bam'
	output: OUTPUT + '{project}/sample_output/{sample}/{sample}_MTB.pos.txt'
	shell:
		"""
		samtools view -F 260 -q 30 {input} | \
			awk '{{ print $4 }}' | \
			sort | uniq -c > {output};
		"""	
