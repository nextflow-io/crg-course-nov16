# Nextflow demo scripts

A set of scripts for Nextflow demo purpose. 


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

## Docker hands-on 

Create a Docker image containing Samtools and Bowtie2. 

Then use it to create a genome index file. 


#### Step 1 

Create an empty working directory eg. `~/docker-demo` and change to it: 

```
mkdir ~/docker-demo 
cd ~/docker-demo 
```

Warning: the Docker build process automatically copies all file in launching working directory 
to the Docker daemon in order to create an image. Thus it's imporant *always* to work in a directory 
containing only the files you want to copy in your Docker image. Alternatively you can use 
the `.dockerignore` file to select the path to exclude from the build. 

#### Step 2 

Use your favourite editor eg. `vim` to create a file named `Dockerfile` and copy the following 
content: 

```

FROM debian:jessie 

MAINTAINER <your name>

RUN apt-get update --fix-missing && \
  apt-get install -q -y python wget unzip samtools
```

When done save the file. 

#### Step 4 

Build the Docker image with the following command: 

```
docker build -t my-image .
```

When it completes, verify that the image has been created listing all available images: 

```
docker images
```

#### Step 5 

Add the Bowtie package to the Docker image adding to the `Dockerfile` the following snippet: 

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

#### Step 6 

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

#### Step 7

Create an genome index file by running bowtie in the container. 

It's important to understand that the container runs a complete separate file system and 
it hasn't access to the hosting file system by default. 

You will need to use the `--volume` command line option to mount the host file system 
in the container. 

Change in the project root directory with the following command: 

```
cd ~/crg-course-nov16
```

The run bowtie in the container with the following command: 

```
docker run --volume $PWD:$PWD --workdir $PWD my-image bowtie2-build data/ggal/genome.fa genome.index
```

#### Step 8 [bonus]

Publish your container in the Docker Hub to share it with other people. 

Create an account in the https://hub.docker.com web site. Then from your shell terminal run 
the following command, entering the user name and password you specified registering in the Hub: 

```
docker login 
``` 

Tag the image with your Docker username account 

```
docker tag my-image <user-name>/my-image 
```

Finally push to the Docker Hub 

```
docker push <user-name>/my-image 
```

After that anyone will be able to download with the command: 

```
docker pull <user-name>/my-image 
```






 
