# Nextflow demo scripts

A set of scripts for Nextflow tutorial purpose. 


## Prerequisite

* Java 7 or 8 
* Docker engine 1.10.x (or higher) 

## Installation 

Install Nextflow by using the following command: 

```
curl -fsSL get.nextflow.io | bash
```
    
Then copy the script `nextflow` in a directory on your `PATH`. 

   
Finally, clone this repository with the following command: 

```
git clone https://github.com/nextflow-io/crg-course-nov16.git && cd crg-course-nov16
```

## Nextflow hands-on 

During this tutorial you will implement a proof of concept of a RNA-Seq pipeline which: 

1. index a genome file
2. map read pairs against the genome
3. perform quantification

### Step 1 

The script `rna-ex1.nf` defines the pipeline input parameters. Run it by using the 
following command: 

```
nextflow run rna-ex1.nf
```

Try to specify a different input parameter, for example: 

```
nextflow run rna-ex1.nf --genome this/and/that
```

### Step 2 

The second example add the `buildIndex` process. It takes the genome file as 
input and create the genome index by using the `bowtie-build` tool. 

Try to run it by using the command: 

```
nextflow run rna-ex2.nf
```

The execution will fail because Bowtie is not installed in the test environment. 

Add the command line option `-with-docker` to launch the execution through a Docker container
as shown below: 

```
nextflow run rna-ex2.nf
```

This it works because it uses the Docker container `nextflow/rnatoy:1.3` defined in the 
`nextflow.config` file. 

In order to avoid to add the option `-with-docker` add the following line in the `nextflow.config` file: 

```
docker.enabled = true
```

### Step 3 

This step shows how to match *read* files into pairs, so that can be mapped by using *TopHat*. 

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

### Step 4 

The script `rna-ex4.nf` add the `mapping` process. Note how it declares three inputs: 
the genome fasta file, the genome index file produced by the `buildIndex` process and 
the read pairs. Also note as the last input is defined as a `set` ie. composed by 
different elements: the pair ID, the first read file and the second read file. 

Execute it by using the following command: 

```
nextflow run rna-ex4.nf -resume
```

The `-resume` option skip the execution of any step that has been processed in a previous 
execution. 

Try to execute it with more read files as shown below: 

```
nextflow run rna-ex4.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq'
```


### Step 5 

This step adds the quantification step to the example script. It takes the 
annotation file and the *bam* files produces by *TopHat* and outputs the transcripts 
*gtf* files. 

You can run it by using the following command: 

```
nextflow run rna-ex5.nf -resume --reads 'data/ggal/reads/*_{1,2}.fq' 
```

### Step 6 

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

You will find the transcripts produces by the pipeline in the `my_transcripts` folder.


### Step 7

This step shows 

### Step 8 (bonus)  

Here you will lean how to publish your pipeline on GitHub. 

Create a new empty project folder eg. `$HOME/rnaseq-demo`. 

Copy in that folder the following files: 

```
cp $HOME/crg-course-nov16/rna-ex6.nf $HOME/rnaseq-demo/main.nf
cp $HOME/crg-course-nov16/nextflow.config $HOME/rnaseq-demo/
cp -r $HOME/crg-course-nov16/data $HOME/rnaseq-demo/
```

Create a new project on GitHub to host your pipeline, and follow 
the instraction provided by it to publish the `$HOME/rnaseq-demo/`
in that repository. 

When done, you will be able to run your pipeline by using the following 
command: 

```
nextflow run <your-github-user-name>/rnaseq-demo
```


## Docker hands-on 

Create a Docker image containing Samtools and Bowtie2. 

Then use it to create a genome index file. 


### Step 1 

Create an empty working directory eg. `~/docker-tutorial` and change to it: 

```
mkdir ~/docker-tutorial 
cd ~/docker-tutorial 
```

Warning: the Docker build process automatically copies all files that are located in the 
launching working directory to the Docker daemon in order to create an image. This can take 
a lot of time for big/many files. For this reason it's imporant *always* to work in a directory 
containing only the files you really need to include in your Docker image. Alternatively you can use 
the `.dockerignore` file to select the path to exclude from the build. 

### Step 2 

Use your favourite editor eg. `vim` to create a file named `Dockerfile` and copy the following 
content: 

```
FROM debian:jessie 

MAINTAINER <your name>

RUN apt-get update --fix-missing && \
  apt-get install -q -y python wget unzip samtools
```

When done save the file. 

### Step 4 

Build the Docker image with the following command: 

```
docker build -t my-image .
```

When it completes, verify that the image has been created listing all available images: 

```
docker images
```

### Step 5 

Add the Bowtie package to the Docker image by adding to the `Dockerfile` the following snippet: 

```
RUN wget -q -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.7/ /opt/bowtie && \
  rm bowtie.zip 

ENV PATH $PATH:/opt/bowtie2-2.2.7/
```

Save the file and build again the image with the same command as before: 

```
docker build -t my-image .
```

You will notice that create a new Docker image with the same *but* with a different image ID. 

### Step 6 

Check that everything is fine running the bowtie help in the container: 

```
docker run my-image bowtie2 --version
```

You can even launch a container in an interactive mode by using the following command: 

```
docker run -it my-image bash
```

Once launched the container you wil noticed that's running as `root`. Use the usual command 
to navigate in the file system. 

To exit from the container, stop the BASH session with the `exit` command.

### Step 7

Create an genome index file by running bowtie in the container. 

The run Bowtie in the container with the following command: 

```
docker run my-image bowtie2-build ~/crg-course-nov16/data/ggal/genome.fa genome.index
```

The above command it fails because Bowtie cannot access the input file.

This is happening because the container runs in a complete separate file system and 
it hasn't access to the hosting file system by default. 

You will need to use the `--volume` command line option to mount the input file(s) eg. 

```
docker run --volume ~/crg-course-nov16/data/ggal/genome.fa:/genome.fa my-image bowtie2-build /genome.fa genome.index
```

An easier way is to mount the parent directory to an equivalent one in the container, 
this allows you to use the same path when running in the container eg. 

```
docker run --volume $HOME:$HOME --workdir $PWD my-image bowtie2-build ~/crg-course-nov16/data/ggal/genome.fa genome.index
```

### Step 8 (bonus)

Publish your container in the Docker Hub to share it with other people. 

Create an account in the https://hub.docker.com web site. Then from your shell terminal run 
the following command, entering the user name and password you specified registering in the Hub: 

```
docker login 
``` 

Tag the image with your Docker username account: 

```
docker tag my-image <user-name>/my-image 
```

Finally push it to the Docker Hub:

```
docker push <user-name>/my-image 
```

After that anyone will be able to download with the command: 

```
docker pull <user-name>/my-image 
```


## Deploy a NF pipeline in the CRG cluster 

Nextflow supports different execution platforms. This means that your script 
can be executed in a single computer or a cluster of by simply providing a configuration
file that specify what resource scheduler you want to use. 

For the sake of this tutorial you will run the [RNA-Toy](https://github.com/nextflow-io/rnatoy) 
pipeline in the CRG cluster. 

Log-in the CRG cluster by using the following cluster: 

```
ssh <sit-xx>@ant-login.linux.crg.es
```

* Replace the `<sit-xx>` string with the user name that you have been assigned. 

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

Then launch the execution of the pipeline with the following command: 

```
nextflow run rnatoy
```

When completed you will find the pipeline output in the `results` folder.

### Error fail over 



### Automatic resource management 

 

## Deploy a NF pipeline in the AWS cloud (bonus)

Nextflow allows any pipeline to the executed in the AWS cloud. 


[![asciicast](https://asciinema.org/a/9vupd4d72ivaz6h56pajjjkop.png)](https://asciinema.org/a/9vupd4d72ivaz6h56pajjjkop)


## Assignment 

Create a two steps pipeline that given any number of protein sequence FASTA files creates 
a phylogenetic tree for each or them. Bonus: use a Docker to isolate and deploy 
the binary dependencies.  

##### Tip 

Use [Clustalw2](http://www.clustal.org/clustal2/) to align the protein sequences. Example 
command line: 

    clustalw2 -infile=sample.fa -output=phylip -outfile=aln.phy
    
Use [RAxML](https://github.com/stamatak/standard-RAxML) to create the phylogenetic tree. 
Example command line: 

    raxmlHPC -f d -j -p 9 -T 2 -m PROTGAMMALG -s aln.phy -n aln       

Use the input protein sequence FASTA files in the following folder: 

    $HOME/crg-course-nov16/
    
    
[The code will be published after the of the course]    

