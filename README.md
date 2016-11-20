# Nextflow + Docker tutorial 

This repository contains the tutorial material for the *Parallel distributed computational workflows
with Nextflow and Docker containers* course. 

## Prerequisite

* Java 7 or 8 
* Docker engine 1.10.x (or higher) 

## Installation 

Install Nextflow by using the following command: 

```
curl -fsSL get.nextflow.io | bash
```
    
The above snippet creates the `nextflow` launcher in the current directory. 
Complete the installation moving it in a directory on your `PATH` eg: 

```
mv nextflow $HOME/bin
``` 
   
Finally, clone this repository with the following command: 

```
git clone https://github.com/nextflow-io/crg-course-nov16.git && cd crg-course-nov16
```

## Nextflow hands-on 

During this tutorial you will implement a proof of concept of a RNA-Seq pipeline which: 

1. Index a genome file.
2. Map read pairs against the genome.
3. Perform quantification.

### Step 1 - Command line parameters

The script `rna-ex1.nf` defines the pipeline input parameters. Run it by using the 
following command: 

```
nextflow run rna-ex1.nf
```

Try to specify a different input parameter, for example: 

```
nextflow run rna-ex1.nf --genome this/and/that
```

### Step 2 - Build genome index

The second example adds the `buildIndex` process. It takes the genome file as 
input and creates the genome index by using the `bowtie-build` tool. 

Try to run it by using the command: 

```
nextflow run rna-ex2.nf
```

The execution will fail because Bowtie is not installed in the test environment. 

Add the command line option `-with-docker` to launch the execution through a Docker container
as shown below: 

```
nextflow run rna-ex2.nf -with-docker
```

This time it works because it uses the Docker container `nextflow/rnatoy:1.3` defined in the 
`nextflow.config` file. 

In order to avoid to add the option `-with-docker` add the following line in the `nextflow.config` file: 

```
docker.enabled = true
```

### Step 3 - Collect read files by pairs

This step shows how to match *read* files into pairs, so thay can be mapped by *TopHat*. 

Edit the script `rna-ex3.nf` and add the following statement as the last line: 

```
read_pairs.println()
```

Save it and execute it with the following command: 

```
nextflow run rna-ex3.nf
```

Try it again specifying different read files by using a glob pattern:

```
nextflow run rna-ex3.nf --reads 'data/ggal/reads/*_{1,2}.fq'
```

It shows how read files matching the pattern specified are grouped in pairs having 
the same prefix.


### Step 4 - Map sequence reads 

The script `rna-ex4.nf` adds the `mapping` process. Note how it declares three inputs: 
the genome fasta file, the genome index file produced by the `buildIndex` process and 
the read pairs. Also note as the last input is defined as a `set` ie. it's composed by 
different elements: the pair ID, the first read file and the second read file. 

Execute it by using the following command: 

```
nextflow run rna-ex4.nf -resume
```

The `-resume` option skips the execution of any step that has been processed in a previous 
execution. 

Try to execute it with more read files as shown below: 

```
nextflow run rna-ex4.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```


### Step 5 - Perform reads quantification 

This step adds the quantification step to the example script. It takes the 
annotation file and the *bam* files produced by *TopHat* and outputs the transcripts 
*gtf* files. 

You can run it by using the following command: 

```
nextflow run rna-ex5.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq' 
```

### Step 6 - Define the pipeline output

This step shows how produce the pipeline output to a folder of your choice by using the 
`publishDir` directive. 

Run the example by using the following command: 

```
nextflow run rna-ex6.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq' 
```

Then you will find the quantification files in the folder `results`. 

Modify the `rna-ex6.nf` script by adding the following line at the beginning of the file: 

```
params.outdir = 'results'
```

Then, look for the `publishDir` directive in the `makeTranscript` process, and 
replace the `'results'` string with the `params.outdir` parameter. 

Finally run it again with the following command: 

```
nextflow run rna-ex6.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq' --outdir my_transcripts
```

You will find the transcripts produced by the pipeline in the `my_transcripts` folder.


### Step 7 - Handle completion event

This step shows how to execute an action when the pipeline completes the execution. 

Note that Nextflow processes define the execution of *asynchronous* tasks i.e. they are not 
executed one after another as they are written in the pipeline script as it would happen in a 
common *iperative* programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message 
when the script completes. 

Try to run it by using the following command: 

```
nextflow run rna-ex7.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
``` 
 
### Step 8 - Manage custom scripts

Real world pipelines use a lot of custom user scripts (BASH, R, Python, etc). Nextflow 
allows you to use and manage all these scripts in consistent manner. Simply put them 
in a directory named `bin` in the pipeline project root. They will be automatically added 
to the pipeline execution `PATH`. 

For example, create a file named `quantify.sh` with the following content: 

```
#!/bin/bash 
set -e 
set -u

annot=${1}
bam_file=${2}
pair_id=${3}

cufflinks --no-update-check -q -G $annot ${bam_file}
mv transcripts.gtf transcript_${pair_id}.gtf
```

Save it, grant the execute permission and move it under the `bin` directory as shown below: 

```
chmod +x quantify.sh
mkdir -p bin 
mv quantify.sh bin
```

Then, open the `rna-ex7.nf` file and replace the `makeTranscript` process with 
the following code: 

```
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'  
       
    input:
    file annot from annotation_file 
    set pair_id, file(bam_file) from bam
     
    output:
    set pair_id, file('transcript_*.gtf') into transcripts
 
    """
    quantify.sh $annot $bam_file $pair_id
    """
}

```

For the sake of simplicity of this example, the *cpus* parameter it's ignored. 

Run it as before: 

```
nextflow run rna-ex7.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```


### Step 9 - Publish to GitHub (bonus)  

Here you will lean how to publish your pipeline on [GitHub](https://github.com) and share 
it with other people and allowing you to track all the project 
dependencies and changes with ease. 
 
Create a new empty project folder eg. 

```
mkdir $HOME/rnaseq-demo
cd $HOME/rnaseq-demo
``` 

Setup your `git` credentials: 

```
git config --config user.name "your name"
git config --config user.email your@email.com 
```

Create a new project on GitHub to host your pipeline, and follow 
the instruction provided by it to publish the project in the project in 
the folder `$HOME/rnaseq-demo/` in that repository. 

Note: make sure to use the same email address you have defined in your 
`git` configuration setup with the previous commands.  

Finally, copy the pipeline files and data and upload to the GitHub repository 

```
cp $HOME/crg-course-nov16/rna-ex6.nf $HOME/rnaseq-demo/main.nf
cp $HOME/crg-course-nov16/nextflow.config $HOME/rnaseq-demo/
cp -r $HOME/crg-course-nov16/bin $HOME/rnaseq-demo/
cp -r $HOME/crg-course-nov16/data $HOME/rnaseq-demo/

git add bin/ data/ main.nf nextflow.config 
git commit -m 'Added pipeline files'
git push 
```

When done, you will be able to run your pipeline by using the following 
command: 

```
nextflow run <your-github-user-name>/rnaseq-demo
```


### Manage revisions (bonus)

Git and GitHub are tools specifically designed to track project changes and versions. 
You can use Git tags, branches or commit IDs to maintain an history revision of your 
pipeline projects.

Nextflow integrates these tools making possible to run any revision of your pipeline 
by simply specifying it on the run command line by using the `-revision` option, as shown 
below: 

```
nextflow run <project name> -r <revision name>
```

The list of available revision can be list by using the following command: 

```
nextflow info <project-name>
```


## Docker hands-on 

Get practice with basic Docker commands to pull, run and build your own containers.
 
A container is a ready-to-run Linux environment which can be executed in an isolated 
manner from the hosting system. It has own copy of the file system, processes space,
memory management, etc. 
 
Containers are a Linux feature known as *Control Groups* or [Ccgroups](https://en.wikipedia.org/wiki/Cgroups)
introduced with kernel 2.6. 

Docker adds to this concept an handy management tool to build, run and share container images. 

These images can be uploaded and published in a centralised repository know as 
[Docker Hub](https://hub.docker.com), or hosted by other parties like for example [Quay](https://quay.io).


### Step 1 - Run a container 

Run a container is easy as using the following command: 

```
docker run <container-name> 
```

For example: 

```
docker run hello-world  
```

### Step 2 - Pull a container 

The pull command allows you to download a Docker image without running it. For example: 

```
docker pull debian:wheezy 
```

The above command download a Debian Linux image.


### Step 3 - Run a container in interactive mode 

Launching a BASH shell in the container allows you to operate in an interactive mode 
in the containerised operating system. For example: 

```
docker run -it debian:wheezy bash 
``` 

Once launched the container you wil noticed that's running as root (!). 
Use the usual commands to navigate in the file system.

To exit from the container, stop the BASH session with the exit command.

### Step 4 - Your first Dockerfile

Docker images are created by using a so called `Dockerfile` i.e. a simple text file 
containing a list of commands to be executed to assemble and configure the image
with the software packages required.    

In this step you will create a Docker image containing the Samtools and Bowtie2 tools.

In order to build a Docker image, start creating an empty directory eg. 
`~/docker-tutorial` and change to it: 

```
mkdir ~/docker-tutorial 
cd ~/docker-tutorial 
```

Warning: the Docker build process automatically copies all files that are located in the 
current directory to the Docker daemon in order to create the image. This can take 
a lot of time when big/many files exist. For this reason it's important to *always* work in 
a directory containing only the files you really need to include in your Docker image. 
Alternatively you can use the `.dockerignore` file to select the path to exclude from the build. 

Then use your favourite editor eg. `vim` to create a file named `Dockerfile` and copy the 
following content: 

```
FROM debian:wheezy 

MAINTAINER <your name>

RUN apt-get update --fix-missing && \
  apt-get install -q -y python wget unzip samtools
```

When done save the file. 


### Step 5 - Build the image  

Build the Docker image by using the following command: 

```
docker build -t my-image .
```

Note: don't miss the dot in the above command. When it completes, verify that the image 
has been created listing all available images: 

```
docker images
```

### Step 6 - Add a software package to the image

Add the Bowtie package to the Docker image by adding to the `Dockerfile` the following snippet: 

```
RUN wget --no-check-certificate -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.7/ /opt/bowtie && \
  rm bowtie.zip 

ENV PATH $PATH:/opt/bowtie2-2.2.7/
```

Save the file and build again the image with the same command as before: 

```
docker build -t my-image .
```

You will notice that it creates a new Docker image with the same name *but* with a 
different image ID. 

### Step 7 - Run Bowtie in the container 

Check that everything is fine running Bowtie in the container as shown below: 

```
docker run my-image bowtie2 --version
```

You can even launch a container in an interactive mode by using the following command: 

```
docker run -it my-image bash
```


### Step 8 - File system mounts

Create an genome index file by running Bowtie in the container. 

Try to run Bowtie in the container with the following command: 

```
docker run my-image \
  bowtie2-build ~/crg-course-nov16/data/ggal/genome.fa genome.index
```

The above command fails because Bowtie cannot access the input file.

This happens because the container runs in a complete separate file system and 
it cannot access the hosting file system by default. 

You will need to use the `--volume` command line option to mount the input file(s) eg. 

```
docker run --volume ~/crg-course-nov16/data/ggal/genome.fa:/genome.fa my-image \
  bowtie2-build /genome.fa genome.index
```

An easier way is to mount a parent directory to an identical one in the container, 
this allows you to use the same path when running it in the container eg. 

```
docker run --volume $HOME:$HOME --workdir $PWD my-image \
  bowtie2-build ~/crg-course-nov16/data/ggal/genome.fa genome.index
```

### Step 9 - Upload the container in the Docker Hub (bonus)

Publish your container in the Docker Hub to share it with other people. 

Create an account in the https://hub.docker.com web site. Then from your shell terminal run 
the following command, entering the user name and password you specified registering in the Hub: 

```
docker login 
``` 

Tag the image with your Docker user name account: 

```
docker tag my-image <user-name>/my-image 
```

Finally push it to the Docker Hub:

```
docker push <user-name>/my-image 
```

After that anyone will be able to download it by using the command: 

```
docker pull <user-name>/my-image 
```


## Deploy a NF pipeline in the CRG cluster 

Nextflow supports different execution platforms. This means that your script 
can be executed in a single computer, a cluster or a cloud by simply providing a configuration
file that specify what computational platform you want to use. 

For the sake of this tutorial you will run the [RNA-Toy](https://github.com/nextflow-io/rnatoy) 
pipeline in the CRG cluster. 

Log-in the CRG cluster by using the following cluster: 

```
ssh <sitXX>@ant-login.linux.crg.es
```

* Replace the `<sitXX>` string with the user name that you have been assigned. 

Create a project directory eg. `rnatoy` and create a file named `nextflow.config` with
the following content: 

```
process.executor = 'crg' 
process.queue = 'course'
process.scratch = true
process.time = '1h'
process.memory = '1G'
docker.enabled = true
```

Then launch the execution of the pipeline by using the following command: 

```
nextflow run rnatoy
```

When completed you will find the pipeline output in the `results` folder.


### Run the pipeline against a real dataset 

Create a new folder to run the pipeline against the mouse genome dataset eg: 

```
mkdir -p $HOME/mouse-run
cd $HOME/mouse-run
```

Then create the `nextflow.config` file with the following content: 

```
params.reads = "/software/rg/rnaseq/data/*_{1,2}.fastq.gz"
params.annot = "/software/rg/rnaseq/refs/mm65.long.ok.sorted.gtf"
params.genome = '/users/cn/ptommaso/nf-course/projects/mouse_genome_mm9_chr1.fa'

process.executor = 'crg' 
process.queue = 'short-sl7'
process.scratch = true
process.time = '1h'
process.memory = '8G'
process.cpus = 4 
process.$buildIndex.cpus = 8 

docker.enabled = true
trace.enabled = true
```

When done, launch the execution by using this command: 

```
nextflow run rnatoy -bg > log
```

The `-bg` will launch NF in the background, to check the execution status you can 
follow the `log` as shown below: 

```
tail -f log
```
 

### Automatic errors fail over 

When running large scale pipelines launching thousands of jobs on many 
different computing nodes errors are not a remote event. 

Nextflow allows failing tasks to be automatically re-executed, in this way it's possible 
to address temporary failures such as failing hardware or network hiccups. In order to enable 
automatic jobs re-execution add the following setting in the `nextflow.config` file: 

```
process.errorStrategy = 'retry'
```

A more common source of errors in computational pipeline are peaks in computing resources, 
allocated by a jobs exceeding the original resource request. In this context automatically 
re-executing the failed task is useless because it would simply replicate the same error condition. 

A common solution consists of increasing the resource request for the needs of the most consuming job, 
even though this will result in a suboptimal allocation of most of the jobs that are less resource hungry.

Nextlow allows resources to be defined in a dynamic manner. In this way it is possible to 
increase the memory request when rescheduling a failing task execution. For example: 

```
process.memory = { 1.GB * task.attempt }
process.errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
```

By using the above settings pipeline a task will initially request one GB of memory. 
In case of an error it will be rescheduled requesting 2 GB and so on, until it is executed 
successfully or the limit of times a task can be retried is reached, forcing the termination 
of the pipeline.


## Deploy a NF pipeline in the AWS cloud (bonus)

Nextflow pipelines can be seamlessly executed in the Amazon cloud. All you need is an AWS 
user account a base Amazon VM image (AMI) that will be used to setup the computing cluster 
in the cloud. 

The following screen cast shows how to configure, setup the cluster and launch the pipeline 
execution in the AWS cloud in a few commands:  

[![asciicast](https://asciinema.org/a/9vupd4d72ivaz6h56pajjjkop.png)](https://asciinema.org/a/9vupd4d72ivaz6h56pajjjkop)


## Assignment 

Create a two steps pipeline that given any number of protein sequence FASTA files creates 
a phylogenetic tree for each or them. Bonus: use a Docker container to isolate and deploy 
the binary dependencies.  

#### Tip 

Use [Clustalw2](http://www.clustal.org/clustal2/) to align the protein sequences. Example 
command line: 

    clustalw2 -infile=sample.fa -output=phylip -outfile=aln.phy
    
Use [RAxML](https://github.com/stamatak/standard-RAxML) to create the phylogenetic tree. 
Example command line: 

    raxmlHPC -f d -j -p 9 -T 2 -m PROTGAMMALG -s aln.phy -n aln       

Use the input protein sequence FASTA files in the following folder: 

    $HOME/crg-course-nov16/data/prot
    
    
[The code will be published after the of the course]    

