# Prep Material - Nextflow and Docker Course
####*November 22, 2016*

Time: 20-25 min

*Developed by Paolo Di Tommaso and Evan Floden*

####The following short exercise should be performed in advance attending of the course.

####*If you have any problems, the first places to visit is the [documentation](https://www.nextflow.io/docs/latest/index.html) or our [Gitter chat](https://gitter.im/nextflow-io/nextflow) where you will find helpful support from the Nextflow community.*

### What is Nextflow?

Nextflow is a *workflow manager*. It has been developed specifically to ease the creation and execution of bioinformatics pipelines.  The benefits of having you pipeline in Nextflow include: 
* Built in GitHub support.
* Automated cluster execution.
* Docker support to eliminate the need for dependencies.
* Portability so you can run your pipeline anywhere laptop, cluster or cloud.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow. 


### What is Docker?

Docker is a tool to package and run applications in isolated environments called containers. In the scientific context, Docker and container technology in general enable us to:
* Simplify the setup of complicated software, libraries and modules.
* Make our results, analysis and graphs 100% reproducible.
* Share our work with the world.


### Why use Nextflow and Docker?

The following narrative will get up and going with Nextflow and Docker. You will need a Unix computer (Linux or Mac) and a little experience using command line applications. 

Far off in the land of *Nonreplicatous*, Marta has just concluded her weekly meeting with her PI Jose. 
Jose has informed her that their paper, so desperately needed to complete her PhD, has been rejected. However, the editor has given them a lifeline if they can submit revisions. Jose explains the new analysis should be easy. A past lab member has previously done a very similar analysis and all the scripts are in his folder. Marta must:

* find the folder
* download the dependencies
* install the software
* modify the scripts
* run the analysis
* upload the results

Finish these tasks and the PhD is finally hers... 

Seven days later and Marta trudges into her weekly meeting the PI. Depressed, she has been stuck for three days trying to install the software. It is incompatible with the new operating system installed on the cluster and she cannot find an old piece of softwar. In that moment, Jose remembers the brilliant summer student who had taken over Mario's analysis. He was very lazy, he hated submitting jobs to the cluster manually and wanted to automate as much as possible.

Jose explained that in his final presentation, the summer student had told them about Nextflow and Docker. He forwards the presentation to Marta in a last ditch attempt to save her PhD. Later that evening and after only two beers at the beer session, the not-so tipsy Marta sits on her sofa and opens the presentation. To her surprise, the presentation contains step-by-step instructions for anyone with a Mac or Linux. She follows along on her laptop:

* Install [Docker](https://www.docker.com/)

    MacOS:
        https://download.docker.com/mac/stable/Docker.dmg
    
    Linux: 
        https://docs.docker.com/engine/installation/
    
Then from the command line/terminal she then installs [Nextflow](http://nextflow.io/index.html#GetStarted)

* Install Nextflow:
    
    ```
    curl -fsSL get.nextflow.io | bash 
    ```

With Nextflow and Docker installed she is then able to run the complete analysis pipeline with two commands.

* Pull the docker container:
   
    ```
    docker pull cbcrg/kallisto-nf
    ```

* Run the analysis:
   
    ```
    nextflow run cbcrg/kallisto-nf -with-docker
    ```
    
Marta arrives in to work the next morning, confident she will obtain the PhD and intrigued at the magic behind the pipeline.
    
## What makes up a Nextflow workflow?

Nextflow introduces two concepts: Processes and Channels. 

**Processes** consist of a defined inputs, outputs and scripts section. Consider our real world kallisto-nf example [here](https://github.com/cbcrg/kallisto-nf/blob/master/kallisto.nf#L76-L90) where the process index takes the transcritome fasta file and generates an index.

**Channels** contain the data and connect processes together.

Both process and channels are defined in the `.nf` script which forms the heart of a pipeline.


## What about Docker?

In Nextflow, each task is sandboxed, effectively executed in an isolated way. By using Docker, each task can be run within its own container instance. Think of it as if every task is executed in its own operating system. The docker container can be defined the in the nextflow config file (`nextflow.config`). 

Consider the kallisto-nf example [here](https://github.com/cbcrg/kallisto-nf/blob/master/nextflow.config#L10). This tells Nextflow to use this Docker container for the execution of each process.


## Summary

The previous introduction provides the basis for how a Nextflow pipeline functions and introduces Docker containers. In class we will see how to create our own workflows and execute them on the cluster. We will also explore how to find or build our own containers and link them to a pipeline.
