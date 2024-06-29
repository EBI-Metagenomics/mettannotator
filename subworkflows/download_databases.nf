/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD DATABASES AUX WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { AMRFINDER_PLUS_GETDB     } from '../modules/local/amrfinder_plus_getdb'
include { ANTISMASH_GETDB          } from '../modules/local/antismash_getdb'
include { DBCAN_GETDB              } from '../modules/local/dbcan_getdb'
include { DEFENSE_FINDER_GETDB     } from '../modules/local/defense_finder_getdb'
include { EGGNOG_MAPPER_GETDB      } from '../modules/local/eggnog_getdb'
include { INTEPROSCAN_GETDB        } from '../modules/local/interproscan_getdb'
include { INTEPRO_ENTRY_LIST_GETDB } from '../modules/local/interpro_list_getdb'
include { RFAM_GETMODELS           } from '../modules/local/rfam_getmodels'
include { BAKTA_GETDB              } from '../modules/local/bakta_getdb'


workflow DOWNLOAD_DATABASES {

    main:

        amrfinder_plus_db = channel.empty()
        antismash_db = channel.empty()
        defense_finder_db = channel.empty()
        dbcan_db = channel.empty()
        interproscan_db = channel.empty()
        interpro_entry_list = channel.empty()
        eggnog_db = channel.empty()
        bakta_db = channel.empty()

        amrfinder_plus_dir = file("$params.dbs/amrfinder/")
        antismash_dir = file("$params.dbs/antismash")
        defense_finder_dir = file("$params.dbs/defense_finder/")
        dbcan_dir = file("$params.dbs/dbcan/")
        interproscan_dir = file("$params.dbs/interproscan")
        interpro_entry_list_dir = file("$params.dbs/interpro_entry_list/")
        eggnog_data_dir = file("$params.dbs/eggnog")
        rfam_ncrna_models = file("$params.dbs/rfam_models/rfam_ncrna_cms")
        bakta_dir = file("$params.dbs/bakta")

        if (amrfinder_plus_dir.exists()) {
            amrfinder_plus_db = tuple(
                amrfinder_plus_dir,
                file("${amrfinder_plus_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
            log.info("AMRFinder plus database exists, or at least the expected folder.")
        } else {
            AMRFINDER_PLUS_GETDB()
            amrfinder_plus_db = AMRFINDER_PLUS_GETDB.out.amrfinder_plus_db.first()
        }

        if (antismash_dir.exists()) {
            antismash_db = tuple(
                antismash_dir,
                file("${antismash_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
            log.info("AntiSMASH database exists, or at least the expected folder.")
        } else {
            ANTISMASH_GETDB()
            antismash_db = ANTISMASH_GETDB.out.antismash_db.first()
        }


        if (defense_finder_dir.exists()) {
            log.info("Defense Finder models exists, or at least the expected folder.")
            defense_finder_dir_db = tuple(
                defense_finder_dir,
                file("${defense_finder_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
        } else {
            DEFENSE_FINDER_GETDB()
            defense_finder_db = DEFENSE_FINDER_GETDB.out.defense_finder_db.first()
        }

        if (dbcan_dir.exists()) {
            log.info("DBCan database exists, or at least the expected folder.")
            dbcan_db = tuple(
                dbcan_dir,
                file("${dbcan_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
        } else {
            DBCAN_GETDB()
            dbcan_db = DBCAN_GETDB.out.dbcan_db.first()
        }

        if (interproscan_dir.exists()) {
            log.info("InterproScan database exists, or at least the expected folder.")
            interproscan_db = tuple(
                interproscan_dir,
                file("${interproscan_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
        } else {
            INTEPROSCAN_GETDB()
            interproscan_db = INTEPROSCAN_GETDB.out.interproscan_db.first()
        }

        if (interpro_entry_list_dir.exists()) {
            log.info("InterPro entry list file exists, or at least the expected folder.")
            interpro_entry_list = tuple(
                interpro_entry_list_dir,
                file("${interpro_entry_list_dir}/VERSION.txt", checkIfExists: true).text // the DB version
            )
        } else {
            INTEPRO_ENTRY_LIST_GETDB()
            interpro_entry_list = INTEPRO_ENTRY_LIST_GETDB.out.interpro_entry_list.first()
        }

        if (eggnog_data_dir.exists()) {
            log.info("EggNOG mapper database exists, or at least the expected folder.")
            eggnog_db = tuple(
                eggnog_data_dir,
                file("${eggnog_data_dir}/VERSION.txt", checkIfExists: true).text
            )
        } else {
            EGGNOG_MAPPER_GETDB()
            eggnog_db = EGGNOG_MAPPER_GETDB.out.eggnog_db.first()
        }

        if (!rfam_ncrna_models.exists()) {
            RFAM_GETMODELS()
            rfam_ncrna_models = RFAM_GETMODELS.out.rfam_ncrna_cms.first()
        } else {
            log.info("RFam model files exists, or at least the expected folders.")
            rfam_ncrna_models = tuple(
                rfam_ncrna_models,
                file("${rfam_ncrna_models}/VERSION.txt", checkIfExists: true).text
            )
        }

        if (params.bakta) {
            if (bakta_db.exists()) {
                log.info("Bakta database exists, or at least the expected folders.")
                bakta_db = tuple(
                    bakta_dir,
                    file("${bakta_dir}/VERSION.txt", checkIfExists: true).text
                )
            } else {
                BAKTA_GETDB()
                bakta_db = BAKTA_GETDB.out.bakta_db.first()
            }
        }
    emit:
        amrfinder_plus_db = amrfinder_plus_db
        antismash_db = antismash_db
        defense_finder_db = defense_finder_db
        dbcan_db = dbcan_db
        interproscan_db = interproscan_db
        interpro_entry_list = interpro_entry_list
        eggnog_db = eggnog_db
        rfam_ncrna_models = rfam_ncrna_models
        bakta_db = bakta_db

}
