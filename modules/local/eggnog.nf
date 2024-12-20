process EGGNOG_MAPPER {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.eggnog-mapper:v2.1.11'

    input:
    // on mode "annotations" will be ignored, submit an empty path (channel.path("NO_FILE"))
    tuple val(meta), file(fasta), file(annotation_hit_table)
    // on mode "mapper" will be ignored, submit an empty path (channel.path("NO_FILE"))
    val mode // mapper or annotations
    tuple path(eggnog_db_dir), val(db_version)

    output:
    tuple val(meta), path("*annotations*"), emit: annotations, optional: true
    tuple val(meta), path("*orthologs*"), emit: orthologs, optional: true
    path "versions.yml", emit: versions

    script:
    def db_mem_flag = ""
    /*
    The required memory of the executor needs to be greater than 44GB
    to be able to load the eggnog SQLite database into memory.
    Docs: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#a-few-recipes
    */
    if ( task.memory >= 44.GB ) {
        db_mem_flag = "--dbmem"
    }
    if ( mode == "mapper" )
        """
        emapper.py -i ${fasta} \
        --database ${eggnog_db_dir}/eggnog.db \
        --dmnd_db ${eggnog_db_dir}/eggnog_proteins.dmnd \
        --data_dir ${eggnog_db_dir} \
        -m diamond \
        --no_file_comments \
        --cpu ${task.cpus} \
        --no_annot ${db_mem_flag} \
        -o ${meta.prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
            eggnog-mapper database: $db_version
        END_VERSIONS
        """
    else if ( mode == "annotations" )
        """
        emapper.py \
        --data_dir ${eggnog_db_dir} \
        --no_file_comments \
        --cpu ${task.cpus} \
        --tax_scope 'prokaryota_broad' \
        --annotate_hits_table ${annotation_hit_table} ${db_mem_flag} \
        -o ${meta.prefix} \
        --override

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
            eggnog-mapper database: $db_version
        END_VERSIONS
        """
    else
        error "Invalid mode: ${mode}"


    // TODO: there must be a more clever way to create the stubs for this.
    stub:
    if ( mode == "mapper" )
        """
        touch eggnog-output.emapper.seed_orthologs

        echo "#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs" > eggnog-output.emapper.seed_orthologs

        echo "MGYG000000012_00001	948106.AWZT01000053_gene1589	1.1e-63	199.0	COG5654@1|root,COG5654@2|Bacteria,1N6P3@1224|Proteobacteria,2VSGY@28216|Betaproteobacteria,1KFU4@119060|Burkholderiaceae	28216|Betaproteobacteria	S	RES	-	-	-	-	-	-	-	-	-	-	-	-	RES" >> eggnog-output.emapper.seed_orthologs

        echo "MGYG000000001_00001	948106.AWZT01000053_gene1589	1.1e-63	199.0	COG5654@1|root,COG5654@2|Bacteria,1N6P3@1224|Proteobacteria,2VSGY@28216|Betaproteobacteria,1KFU4@119060|Burkholderiaceae	28216|Betaproteobacteria	S	RES	-	-	-	-	-	-	-	-	-	-	-	-	RES" >> eggnog-output.emapper.seed_orthologs

        echo "MGYG000000020_00001	948106.AWZT01000053_gene1589	1.1e-63	199.0	COG5654@1|root,COG5654@2|Bacteria,1N6P3@1224|Proteobacteria,2VSGY@28216|Betaproteobacteria,1KFU4@119060|Burkholderiaceae	28216|Betaproteobacteria	S	RES	-	-	-	-	-	-	-	-	-	-	-	-	RES" >> eggnog-output.emapper.seed_orthologs
        """
    else if ( mode == "annotations" )
        """
        touch eggnog-output.emapper.annotations
        echo "#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs" > eggnog-output.emapper.annotations

        echo "MGYG000000012_00001	59538.XP_005971304.1	7.97e-152	431.0	COG0101@1|root,KOG4393@2759|Eukaryota,39RAQ@33154|Opisthokonta,3BK4Y@33208|Metazoa,3D27W@33213|Bilateria,48A93@7711|Chordata,494G6@7742|Vertebrata,3J2WS@40674|Mammalia 33208|Metazoa	J	synthase-like 1 -	GO:0001522	-	-	-	-	-	-	-	-	-	-	DSPc,Laminin_G_3,PseudoU_synth_1" >> eggnog-output.emapper.annotations

        echo "MGYG000000001_00001	 59538.XP_005971304.1	7.97e-152	431.0	COG0101@1|root,KOG4393@2759|Eukaryota,39RAQ@33154|Opisthokonta,3BK4Y@33208|Metazoa,3D27W@33213|Bilateria,48A93@7711|Chordata,494G6@7742|Vertebrata,3J2WS@40674|Mammalia 33208|Metazoa	J	synthase-like 1 -	GO:0001522	-	-	-	-	-	-	-	-	-	-	DSPc,Laminin_G_3,PseudoU_synth_1" >> eggnog-output.emapper.annotations

        echo "MGYG000000020_00001	59538.XP_005971304.1	7.97e-152	431.0	COG0101@1|root,KOG4393@2759|Eukaryota,39RAQ@33154|Opisthokonta,3BK4Y@33208|Metazoa,3D27W@33213|Bilateria,48A93@7711|Chordata,494G6@7742|Vertebrata,3J2WS@40674|Mammalia 33208|Metazoa	J	synthase-like 1 -	GO:0001522	-	-	-	-	-	-	-	-	-	-	DSPc,Laminin_G_3,PseudoU_synth_1" >> eggnog-output.emapper.annotations
        """
    else
        error "Invalid mode: ${mode}"
}
