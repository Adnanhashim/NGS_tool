#!/bin/bash

# Constants
DIR="/cluster/projects/nn8009k/Adnan"
SCRIPTS="$DIR/000_NGS_tool"

# Function to create directories if they don't exist
create_directory() {
  [ -d "$1" ] || mkdir -p "$1"
}

# Function to check if a directory exists or create it
check_or_create_directory() {
  local dir="$1"
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
    if [ $? -ne 0 ]; then
      echo "Error: Failed to create directory $dir"
      exit 1
    fi
  fi
}


# Function to perform fastqc analysis

run_fastqc() {
    # Check if the required variables are set
        if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
          echo "Some or all of the required variables are not set properly."
          return
        fi

    # Define a log file for recording errors
        local log_file="fastqc_error_log.txt"

    # Create or truncate the log file
        #> "$log_file"
        > "$Output_dir/$log_file"

    # Define a function for error handling
        handle_error() {
          local error_message="$1"
          echo "Error: $error_message"
          echo "Error: $error_message" >> "$Output_dir/$log_file"
        }

  # Define the output directory based on the layout
        if [ "$layout" = "PE" ]; then
            fastqc_dir="$Output_dir/OUTPUT/1_fastqc_PE"
        elif [ "$layout" = "SE" ]; then
            fastqc_dir="$Output_dir/OUTPUT/1_fastqc_SE"
        else
            handle_error "Invalid layout: $layout"
            return
        fi

  # Create the output directory if it doesn't exist
        if [ ! -e "$fastqc_dir" ]; then
          mkdir -p "$fastqc_dir"
          if [ $? -ne 0 ]; then
            handle_error "Failed to create the output directory: $fastqc_dir"
            return
          fi
        fi

  # Loop through input FastQ files and run FastQC
      for fq in "$Input_dir"/*.fastq.gz; do
          if [ ! -e "$fq" ]; then
              handle_error "Input file $fq does not exist."
              continue
          fi

          # Run FastQC on the input file
          echo "sbatch $SCRIPTS/1_fastqc.sh $fq $fastqc_dir"
          sbatch "$SCRIPTS/1_fastqc.sh" "$fq" "$fastqc_dir"

          # Check if sbatch executed successfully
          if [ $? -ne 0 ]; then
              handle_error "Failed to execute sbatch for file: $fq"
          fi
      done

  # Check if there were any errors recorded in the log file
      if [ -s "$log_file" ]; then
         echo "FastQC analysis completed with errors. See $log_file for details."
      else
        echo "FastQC analysis commands submitted successfully."
      fi
}


# Function to perform Trimglore analysis
run_trimglore() {
    # Check if the required variables are set
        if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
          echo "Some or all of the required variables are not set properly."
          return
        fi

    # Define a log file for recording errors
        local log_file="trimglore_error_log.txt"

    # Create or truncate the log file
        > "$Output_dir/$log_file"

    # Define a function for error handling
        handle_error() {
            local error_message="$1"
            echo "Error: $error_message"
            echo "Error: $error_message" >> "$Output_dir/$log_file"
        }

    # Define the output directory based on the layout
        if [ "$layout" = "PE" ]; then
            trimglore_dir="$Output_dir/OUTPUT/2_Trimglore/PE"
            fastqc_dir="$trimglore_dir/fastqc"  # Define the fastqc directory within the layout-specific directory
        elif [ "$layout" = "SE" ]; then
            trimglore_dir="$Output_dir/OUTPUT/2_Trimglore/SE"
            fastqc_dir="$trimglore_dir/fastqc" 
        else
            handle_error "Invalid layout: $layout"
            return
        fi

    # Create the output directory if it doesn't exist
        if [ ! -e "$trimglore_dir" ]; then
          mkdir -p "$trimglore_dir"
            if [ $? -ne 0 ]; then
                handle_error "Failed to create the output directory: $trimglore_dir"
                return
            fi
        fi

    # Create the fastqc directory within the layout-specific directory
        if [ ! -e "$fastqc_dir" ]; then
          mkdir -p "$fastqc_dir"
            if [ $? -ne 0 ]; then
              handle_error "Failed to create the fastqc directory: $fastqc_dir"
              return
            fi
        fi


    # Loop through input FastQ files and run Trimglore
        for fq1 in "$Input_dir"/*_R1*fastq.gz; do
            fq2="${fq1/_R1_/_R2_}"
            out="$(basename "$fq1" _R1_001.fastq.gz)"

            if [ ! -e "$fq1" ] || [ ! -e "$fq2" ]; then
                handle_error "Error: R1 and R2 files do not match for $fq1 and $fq2"
                continue
            fi

    # Run Trimglore on the input files
          echo "sbatch $SCRIPTS/2_trimglore_PE.sh $fq1 $fq2 $trimglore_dir $out"
          sbatch "$SCRIPTS/2_trimglore_PE.sh" "$fq1" "$fq2" "$trimglore_dir" "$out"

    # Check if sbatch executed successfully
          if [ $? -ne 0 ]; then
            handle_error "Failed to execute sbatch for file: $fq1"
          fi
      done

    # Check if there were any errors recorded in the log file
        if [ -s "$log_file" ]; then
            echo "Trimglore analysis completed with errors. See $log_file for details."
        else
            echo "Trimglore analysis commands submitted successfully."
        fi
}


    # Function to perform Alignment analysis
          run_alignment() {
    
    # Check if the required variables are set
        if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
            echo "Some or all of the required variables are not set properly."
            return
        fi

    # Define a log file for recording errors
        local log_file="alignment_error_log.txt"

    # Create or truncate the log file
        > "$Output_dir/$log_file"

    # Define a function for error handling
      handle_error() {
        local error_message="$1"
        echo "Error: $error_message"
        echo "Error: $error_message" >> "$Output_dir/$log_file"
      }

    # Define the output directory based on the layout
        if [ "$layout" = "PE" ]; then
            alignment_dir="$Output_dir/OUTPUT/3_bowtie2/PE"
        elif [ "$layout" = "SE" ]; then
            alignment_dir="$Output_dir/OUTPUT/3_bowtie2/SE"
        else
            handle_error "Invalid layout: $layout"
            return
        fi

    # Create the output directory if it doesn't exist
        if [ ! -e "$alignment_dir" ]; then
          mkdir -p "$alignment_dir"
            if [ $? -ne 0 ]; then
              handle_error "Failed to create the output directory: $alignment_dir"
              return
            fi
        fi
#####################################
  # Loop through input FastQ files and run Bowtie2 alignment
  if [ "$layout" = "PE" ]; then
    for fq1 in "$Input_dir"/*_R1_001_val_1.fq.gz; do
      fq2="${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
      out="$(basename "$fq1" _R1_001_val_1.fq.gz)"

      if [ ! -e "$fq1" ] || [ ! -e "$fq2" ]; then
        handle_error "Error: R1 and R2 files do not match for $fq1 and $fq2"
        continue
      fi

      # Run Bowtie2 alignment on the input files
      echo "sbatch $SCRIPTS/3_bowtie2_PE.sh $fq1 $fq2 $alignment_dir $specie"
      sbatch "$SCRIPTS/3_bowtie2_PE.sh" "$fq1" "$fq2" "$alignment_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  elif [ "$layout" = "SE" ]; then
    for fq1 in "$Input_dir"/*trimmed.fq.gz; do
      out="$(basename "$fq1" trimmed.fq.gz)"

      # Run Bowtie2 alignment on the input files
      echo "sbatch $SCRIPTS/3_bowtie2_SE.sh $fq1 $alignment_dir $specie"
      sbatch "$SCRIPTS/3_bowtie2_SE.sh" "$fq1" "$alignment_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  fi

  # Check if there were any errors recorded in the log file
  if [ -s "$Output_dir/$log_file" ]; then
    echo "Bowtie2 alignment completed with errors. See $Output_dir/$log_file for details."
  else
    echo "Bowtie2 alignment commands submitted successfully."
  fi
}

# Function to perform Hisat2 Alignment analysis

run_hisat2() {
  # Check if the required variables are set
  if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
    echo "Some or all of the required variables are not set properly."
    return
  fi

  # Define a log file for recording errors
  local log_file="hisat2_error_log.txt"

  # Create or truncate the log file
  > "$Output_dir/$log_file"

  # Define a function for error handling
  handle_error() {
    local error_message="$1"
    echo "Error: $error_message"
    echo "Error: $error_message" >> "$Output_dir/$log_file"
  }

  # Define the output directory based on the layout
  if [ "$layout" = "PE" ]; then
    hisat2_dir="$Output_dir/OUTPUT/3_Hisat2/PE"
  elif [ "$layout" = "SE" ]; then
    hisat2_dir="$Output_dir/OUTPUT/3_Hisat2/SE"
  else
    handle_error "Invalid layout: $layout"
    return
  fi

  # Create the output directory if it doesn't exist
  if [ ! -e "$hisat2_dir" ]; then
    mkdir -p "$hisat2_dir"
    if [ $? -ne 0 ]; then
      handle_error "Failed to create the output directory: $hisat2_dir"
      return
    fi
  fi

  # Loop through input FastQ files and run Hisat2
  if [ "$layout" = "PE" ]; then
    for fq1 in "$Input_dir"/*_R1_001_val_1.fq.gz; do
      fq2="${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
      out="$(basename "$fq1" _R1_001_val_1.fq.gz)"

      if [ ! -e "$fq1" ] || [ ! -e "$fq2" ]; then
        handle_error "Error: R1 and R2 files do not match for $fq1 and $fq2"
        continue
      fi

      # Run Hisat2 on the input files
      echo "sbatch $SCRIPTS/hisat2_PE.sh $fq1 $fq2 $hisat2_dir $specie"
      sbatch "$SCRIPTS/hisat2_PE.sh" "$fq1" "$fq2" "$hisat2_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  elif [ "$layout" = "SE" ]; then
    for fq1 in "$Input_dir"/*trimmed.fq.gz; do
      out="$(basename "$fq1" trimmed.fq.gz)"

      # Run Hisat2 on the input files
      echo "sbatch $SCRIPTS/hisat2_SE.sh $fq1 $hisat2_dir $specie"
      sbatch "$SCRIPTS/hisat2_SE.sh" "$fq1" "$hisat2_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  fi

  # Check if there were any errors recorded in the log file
  if [ -s "$Output_dir/$log_file" ]; then
    echo "Hisat2 analysis completed with errors. See $Output_dir/$log_file for details."
  else
    echo "Hisat2 analysis commands submitted successfully."
  fi
}


# Function to perform DNA meth Alignment analysis
run_bismark() {
  # Check if the required variables are set
  if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
    echo "Some or all of the required variables are not set properly."
    return
  fi

  # Define a log file for recording errors
  local log_file="bismark_error_log.txt"

  # Create or truncate the log file
  > "$Output_dir/$log_file"

  # Define a function for error handling
  handle_error() {
    local error_message="$1"
    echo "Error: $error_message"
    echo "Error: $error_message" >> "$Output_dir/$log_file"
  }

  # Define the output directory based on the layout
  if [ "$layout" = "PE" ]; then
    bismark_dir="$Output_dir/OUTPUT/3_bismark/PE"
  elif [ "$layout" = "SE" ]; then
    bismark_dir="$Output_dir/OUTPUT/3_bismark/SE"
  else
    handle_error "Invalid layout: $layout"
    return
  fi

  # Create the output directory if it doesn't exist
  if [ ! -e "$bismark_dir" ]; then
    mkdir -p "$bismark_dir"
    if [ $? -ne 0 ]; then
      handle_error "Failed to create the output directory: $bismark_dir"
      return
    fi
  fi

  # Loop through input FastQ files and run Bismark
  if [ "$layout" = "PE" ]; then
    for fq1 in "$Input_dir"/*_R1_001_val_1.fq.gz; do
      fq2="${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
      out="$(basename "$fq1" _R1_001_val_1.fq.gz)"

      if [ ! -e "$fq1" ] || [ ! -e "$fq2" ]; then
        handle_error "Error: R1 and R2 files do not match for $fq1 and $fq2"
        continue
      fi

      # Run Bismark on the input files
      echo "sbatch $SCRIPTS/bismark_PE.sh $fq1 $fq2 $bismark_dir $specie"
      sbatch "$SCRIPTS/bismark_PE.sh" "$fq1" "$fq2" "$bismark_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  elif [ "$layout" = "SE" ]; then
    for fq1 in "$Input_dir"/*trimmed.fq.gz; do
      out="$(basename "$fq1" trimmed.fq.gz)"

      # Run Bismark on the input files
      echo "sbatch $SCRIPTS/bismark_SE.sh $fq1 $bismark_dir $specie"
      sbatch "$SCRIPTS/bismark_SE.sh" "$fq1" "$bismark_dir" "$specie"

      # Check if sbatch executed successfully
      if [ $? -ne 0 ]; then
        handle_error "Failed to execute sbatch for file: $fq1"
      fi
    done
  fi

  # Check if there were any errors recorded in the log file
  if [ -s "$Output_dir/$log_file" ]; then
    echo "Bismark analysis completed with errors. See $Output_dir/$log_file for details."
  else
    echo "Bismark analysis commands submitted successfully."
  fi
}


# Function to perform ATAC seq Alignment analysis
run_atac() {
  # Check if the required variables are set
  if [ -z "$Module" ] || [ -z "$layout" ] || [ -z "$Input_dir" ] || [ -z "$Output_dir" ]; then
    echo "Some or all of the required variables are not set properly."
    return
  fi

  # Define a log file for recording errors
  local log_file="atac_error_log.txt"

  # Create or truncate the log file
  > "$Output_dir/$log_file"

  # Define a function for error handling
  handle_error() {
    local error_message="$1"
    echo "Error: $error_message"
    echo "Error: $error_message" >> "$Output_dir/$log_file"
  }

  if [ "$Module" = "ATAC" ]; then
    if [ "$layout" = "PE" ] || [ "$layout" = "SE" ]; then
      if [ ! -e "$Output_dir" ]; then
        mkdir -p "$Output_dir/OUTPUT/3_bowtie2_ATAC/$layout"
      fi

      for fq1 in $(ls "$Input_dir"/*${layout}_R1_001_val_1.fq.gz); do
        fq2="${fq1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
        out=$(basename "$fq1" _R1_001_val_1.fq.gz)

        if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
          handle_error "Error: R1 and R2 files do not match for $out"
          continue
        fi

        # Run Bowtie2 on the input files
        echo "sbatch $SCRIPTS/3_bowtie2_ATAC_${layout}.sh $fq1 $fq2 $Output_dir/OUTPUT/3_bowtie2_ATAC/$layout $specie"
        sbatch "$SCRIPTS/3_bowtie2_ATAC_${layout}.sh" "$fq1" "$fq2" "$Output_dir/OUTPUT/3_bowtie2_ATAC/$layout" "$specie"

        # Check if sbatch executed successfully
        if [ $? -ne 0 ]; then
          handle_error "Failed to execute sbatch for file: $out"
        fi
      done
    else
      handle_error "Invalid layout: $layout. Please use 'PE' or 'SE'."
    fi
  else
    handle_error "Please select the right module (ATAC) and layout (PE or SE)."
  fi

  # Check if there were any errors recorded in the log file
  if [ -s "$Output_dir/$log_file" ]; then
    echo "ATAC analysis completed with errors. See $Output_dir/$log_file for details."
  else
    echo "ATAC analysis commands submitted successfully."
  fi
}
