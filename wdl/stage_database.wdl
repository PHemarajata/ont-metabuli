version 1.0

## Stage a Metabuli database into your Terra workspace bucket — run ONCE.
##
## The main workflow takes the database as a `File` (a .tar.gz in GCS) so that
## Cromwell can localize and call-cache it. This helper downloads a prebuilt
## database straight into the workspace bucket, which avoids pulling ~74 GB
## down to a laptop and pushing it back up again.
##
## After it finishes, copy the `database_tar` output path (gs://...) and use it
## as `ont_metabuli.metabuli_db_tar` from then on.
##
## Databases (https://opendata.mmseqs.org/metabuli) all require Metabuli >=1.2.0:
##   refseq_standard  ~75 GB   archaea, bacteria, virus, plasmid, protozoa, fungi + human
##   hrom             ~22 GB   human oral microbiome + human + virus
##   hrgm2            ~85 GB   human gut microbiome v2 + human + virus
##   gtdb226/gtdb232  ~378 GB  full GTDB + human + virus
##
## NOTE: there is no small viral-only prebuilt database any more; the former
## RefSeq_virus (~4 GB) has been withdrawn from every host.

workflow stage_metabuli_database {
    input {
        String db_name    = "refseq_standard"
        String db_host    = "https://opendata.mmseqs.org/metabuli"
        Int    disk_gb    = 200
        Int    cpu        = 4
        Int    memory_gb  = 8
    }

    call DownloadDatabase {
        input:
            db_name   = db_name,
            db_host   = db_host,
            disk_gb   = disk_gb,
            cpu       = cpu,
            memory_gb = memory_gb
    }

    output {
        File   database_tar = DownloadDatabase.database_tar
        String database_size = DownloadDatabase.size_report
    }
}

task DownloadDatabase {
    input {
        String db_name
        String db_host
        Int    disk_gb
        Int    cpu
        Int    memory_gb
    }

    command <<<
        set -euo pipefail
        URL="~{db_host}/~{db_name}.tar.gz"
        echo "downloading $URL"

        # -C - resumes if the task is retried after a partial transfer
        curl -fSL -C - --retry 10 --retry-delay 10 --retry-connrefused \
            "$URL" -o "~{db_name}.tar.gz"

        # fail loudly rather than emitting an HTML/XML error page as a database
        if [ "$(file -b --mime-type "~{db_name}.tar.gz")" != "application/gzip" ]; then
            echo "ERROR: download is not gzip — the URL may have moved:" >&2
            head -c 300 "~{db_name}.tar.gz" >&2
            exit 1
        fi
        gzip -t "~{db_name}.tar.gz"

        du -h "~{db_name}.tar.gz" | cut -f1 > size.txt
        echo "downloaded $(cat size.txt)"
    >>>

    output {
        File   database_tar = "~{db_name}.tar.gz"
        String size_report  = read_string("size.txt")
    }
    runtime {
        docker: "google/cloud-sdk:slim"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: 0
        maxRetries: 2
    }
}
