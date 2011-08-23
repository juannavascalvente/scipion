#!/usr/bin/env python
import glob
import os
import sys
from operator import itemgetter
from threading import Thread
from protlib_utils import reportError
xmippdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]
sys.path.append(xmippdir+'/lib') # add default search path
sys.path.append(xmippdir+'/protocols') # add default search path
import xmipp
#import launch_job

# --------------------------------------------------------------------
def removeOutliersWithinCluster(Nproc,fnRoot,classBlockRefs,thZscore,thPCAZscore,
    analyzingCore,systemFlavour):
    print "Analyzing within cluster ------------------------------------------"
    fh_mpi=os.open(fnRoot + '_multivariate_analysis.sh', os.O_WRONLY | os.O_TRUNC | os.O_CREAT, 0700)
    PCAoutputFiles=[]
    outputFiles=[]
    retval=[]
    for classBlockRef in classBlockRefs:
        blockName=classBlockRef[0]
        classFile=classBlockRef[1]
        reference=classBlockRef[2]
        classDir,classRoot=os.path.split(classFile.replace(".sel",""))

        PCAoutputFile=classDir+"/PCAscore_cluster_"+blockName+"_"+classRoot+".sel"
        command="xmipp_classify_analyze_cluster"+\
                " -i "+blockName+"@"+classFile+\
                " --ref "+reference+\
                " -o "+PCAoutputFile+\
                " --maxDist "+str(thPCAZscore)+\
                " --quiet"
        if analyzingCore:
            if not os.path.exists(classDir+"/Basis"):
                os.mkdir(classDir+"/Basis")
            alignedFile=classDir+"/Basis/PCAscore_cluster_"+blockName+"_"+classRoot+"_aligned.stk"
            basisFile=classDir+"/Basis/PCAscore_cluster_"+blockName+"_"+classRoot+"_basis.stk"
            retval.append((blockName,alignedFile,basisFile,PCAoutputFile))
            command+=" --produceAligned "+alignedFile+\
                     " --basis "+basisFile
        PCAoutputFiles.append(PCAoutputFile)
        os.write(fh_mpi, command+"\n");

        outputFile=classDir+"/score_cluster_"+blockName+"_"+classRoot
        command="xmipp_image_sort"+\
                " -i "+blockName+"@"+classFile+\
                " -o "+outputFile+\
                " --multivariate -v 0"
        outputFiles.append(outputFile+".sel")
        os.write(fh_mpi, command+"\n");
    
    os.close(fh_mpi)
    os.system(launch_job.buildCommand('xmipp_run',
             '-i '+fnRoot + '_multivariate_analysis.sh',True,Nproc,1,systemFlavour))
    os.remove(fnRoot + '_multivariate_analysis.sh')
    
    # Pickup results
    Zscore={}
    Nscore={}
    for outputFile in outputFiles:
        MD=xmipp.MetaData(outputFile)
        for id in MD:
            image=MD.getValue(xmipp.MDL_IMAGE)
            zscore=MD.getValue(xmipp.MDL_ZSCORE)
            if image in Zscore:
                Zscore[image]=Zscore[image]+zscore
                Nscore[image]=Nscore[image]+1
            else:
                Zscore[image]=zscore
                Nscore[image]=1
        os.remove(outputFile)
    PCAZscore={}
    PCANscore={}
    for outputFile in PCAoutputFiles:
        MD=xmipp.MetaData(outputFile)
        for id in MD:
            if analyzingCore:
                image=MD.getValue(xmipp.MDL_IMAGE_ORIGINAL)
            else:
                image=MD.getValue(xmipp.MDL_IMAGE)
            PCAzscore=MD.getValue(xmipp.MDL_ZSCORE)
            if image in PCAZscore:
                PCAZscore[image]=PCAZscore[image]+PCAzscore
                PCANscore[image]=PCANscore[image]+1
            else:
                PCAZscore[image]=PCAzscore
                PCANscore[image]=1
        if not analyzingCore:
            os.remove(outputFile)
    
    # Compute average Zscore and PCAZscore
    for image in Zscore:
        Zscore[image]=Zscore[image]/Nscore[image]
    for image in PCAZscore:
        PCAZscore[image]=PCAZscore[image]/PCANscore[image]

    # Filter classes
    for classBlockRef in classBlockRefs:
        MDin=xmipp.MetaData()
        blockName=classBlockRef[0]
        classFile=classBlockRef[1]
        MDin.readBlock(classFile,blockName)
        MDout=xmipp.MetaData()
        for id in MDin:
            image=MDin.getValue(xmipp.MDL_IMAGE)
            if Zscore[image]<thZscore and PCAZscore[image]<thPCAZscore:
                MDout.addObject()
                MDout.setValue(xmipp.MDL_IMAGE,image)
                MDout.setValue(xmipp.MDL_ZSCORE,(Zscore[image]+PCAZscore[image])/2)
        MDout.writeBlock(classFile,blockName+"_filtered")

    # Update final classification
    if not analyzingCore:
        classFile=fnRoot+".sel"
        stackFile=classFile.replace(".sel",".stk")
        numberOfCodes=xmipp.SingleImgSize(stackFile)[3]
        for codeNumber in range(numberOfCodes):
            MDin=xmipp.MetaData()
            MDin.readBlock(classFile,"class_%06d"%codeNumber)
            MDout=xmipp.MetaData()
            for id in MDin:
                image=MDin.getValue(xmipp.MDL_IMAGE)
                imagePreprocessed=image.replace('class_aligned','preprocessedImages')
                if imagePreprocessed in Zscore:
                    if Zscore[imagePreprocessed]<thZscore and PCAZscore[imagePreprocessed]<thPCAZscore:
                        MDout.addObject()
                        MDout.setValue(xmipp.MDL_IMAGE,image)
                        MDout.setValue(xmipp.MDL_ZSCORE,
                            (Zscore[imagePreprocessed]+PCAZscore[imagePreprocessed])/2)
            MDout.writeBlock(classFile,"class_%06d"%codeNumber+"_filtered")
    return retval
    
# --------------------------------------------------------------------
class buildCoreThread(Thread):
    def __init__(self,Nproc,myId,fnRoot,level,Nclasses,\
        thGoodClass,Nlevels):
        Thread.__init__(self)
        self.Nproc=Nproc
        self.myId=myId
        self.fnRoot=fnRoot
        self.level=level
        self.Nlevels=Nlevels
        self.Nclasses=Nclasses
        self.thGoodClass=thGoodClass
        self.totalSum=0
        self.goodSum=0
        self.classSum=0
        self.allCores=set()

    def run(self):
        classFile=self.fnRoot+"_level_%02d.sel"%self.level
        for nclass in range(0,self.Nclasses):
            if ((nclass+1) % self.Nproc)==self.myId:
                # Read the selfile for this class and level
                blockName="class_%06d"%nclass
                thisClass=xmipp.MetaData()
                thisClass.readBlock(classFile,blockName+"_filtered")
                self.totalSum=self.totalSum+thisClass.size()

                # For all previous levels
                coocurrence={}
                for levelp in range(0,self.level):
                    classFilep=self.fnRoot+"_level_%02d.sel"%levelp
                    # Get the number of classes in that level
                    Nclasses=xmipp.SingleImgSize(fnRoot+"_level_%02d.stk"%levelp)[3]
                    for nclassp in range(0,Nclasses):
                        blockNamep="class_%06d"%nclassp
                        otherClass=xmipp.MetaData()
                        otherClass.readBlock(classFilep,blockNamep+"_filtered")
                        otherClass.intersection(thisClass,xmipp.MDL_IMAGE)
                        
                        if (otherClass.size()>1):
                            filesRemaining=[]
                            for id in otherClass:
                                filesRemaining.append(otherClass.getValue(xmipp.MDL_IMAGE))
                            Nremaining=len(filesRemaining)
                            for id1 in range(Nremaining):
                                for id2 in range(Nremaining):
                                    if id2>id1:
                                        key=(filesRemaining[id1],filesRemaining[id2])
                                        if key in coocurrence:
                                            coocurrence[key]=coocurrence[key]+1
                                        else:
                                            coocurrence[key]=1

                # Analyze coocurrence
                core=set()
                if len(coocurrence)>0:
                    for key in coocurrence:
                        if coocurrence[key]==self.level:
                            core.add(key[0])
                            core.add(key[1])
                if thisClass.size()>0:
                    if (float(len(core))/thisClass.size())>self.thGoodClass:
                        self.goodSum=self.goodSum+len(core)
                        self.classSum=self.classSum+1
                    else:
                        core=set()
                        
                    MDout=xmipp.MetaData()
                    for image in core:
                        MDout.addObject()
                        MDout.setValue(xmipp.MDL_IMAGE,image)
                    MDout.sort(xmipp.MDL_IMAGE)
                    MDout.writeBlock(classFile,blockName+"_core")
                    if self.level==(self.Nlevels-1):
                        self.allCores=self.allCores.union(core)

def computeCores(Nproc,fnRoot,Nlevels):
    print "Computing cores -------------------------------------------"
    allCoresFinalLevel=set()
    # for level in range(1,Nlevels):
    if Nlevels==1:
        level=0

        # Get the number of classes in that level
        Nclasses=xmipp.SingleImgSize(fnRoot+"_level_%02d.stk"%level)[3]
        print "Analyzing "+str(Nclasses)+" classes at level "+str(level)

        # Analyze this level
        classFile=fnRoot+"_level_%02d.sel"%level
        totalSum=0
        for nclass in range(0,Nclasses):
            # Read the selfile for this class and level
            blockName="class_%06d"%nclass
            thisClass=xmipp.MetaData()
            thisClass.readBlock(classFile,blockName+"_filtered")
            thisClass.writeBlock(classFile,blockName+"_core")
            totalSum=totalSum+thisClass.size()
            for id in thisClass:
                fnImg=thisClass.getValue(xmipp.MDL_IMAGE)
                allCoresFinalLevel.add(fnImg)
        goodSum=totalSum
        classSum=Nclasses

        print "   "+str(goodSum)+" out of "+str(totalSum)+\
            " images are in "+str(classSum)+" good classes"
    else:
        for level in range(Nlevels-1,Nlevels):
            # Get the number of classes in that level
            Nclasses=xmipp.SingleImgSize(fnRoot+"_level_%02d.stk"%level)[3]

            # Analyze this level in threads
            print "Analyzing "+str(Nclasses)+" classes at level "+str(level)
            threadList=[]
            for i in range(0,Nproc):
                currentThread=buildCoreThread(Nproc,i,
                    fnRoot,level,Nclasses,thGoodClass,Nlevels)
                threadList.append(currentThread)
                currentThread.start()

            # Pickup results from threads
            totalSum=0
            goodSum=0
            classSum=0
            for thread in threadList:
                thread.join()
                totalSum=totalSum+thread.totalSum
                goodSum=goodSum+thread.goodSum
                classSum=classSum+thread.classSum
                if level==(Nlevels-1):
                    allCoresFinalLevel=allCoresFinalLevel.union(thread.allCores)

            print "   "+str(goodSum)+" out of "+str(totalSum)+\
                " images are in "+str(classSum)+" good classes"
        
    # Update final classification
    classFile=fnRoot+".sel"
    stackFile=classFile.replace(".sel",".stk")
    numberOfCodes=xmipp.SingleImgSize(stackFile)[3]
    for codeNumber in range(numberOfCodes):
        MDin=xmipp.MetaData()
        MDin.readBlock(classFile,"class_%06d"%codeNumber+"_filtered")
        MDout=xmipp.MetaData()
        for id in MDin:
            image=MDin.getValue(xmipp.MDL_IMAGE)
            imagePreprocessed=image.replace('class_aligned','preprocessedImages')
            if imagePreprocessed in allCoresFinalLevel:
                MDout.addObject()
                MDout.setValue(xmipp.MDL_IMAGE,image)
        MDout.writeBlock(classFile,"class_%06d"%codeNumber+"_core")        

# --------------------------------------------------------------------
def gatherCores(fnRoot,classFile,coreInformation):
    MDbasis=xmipp.MetaData()
    MDclassAverages=xmipp.MetaData()
    MDclassAveragesAux=xmipp.MetaData()
    classIdx={}
    idx=0
    for id in coreInformation:
        block=id[0]
        alignedStack=id[1]
        basisStack=id[2]
        PCAoutput=id[3]
        
        # Put all the basis together
        numberOfImages=xmipp.SingleImgSize(basisStack)[3]
        for i in range(numberOfImages):
            MDbasis.addObject()
            MDbasis.setValue(xmipp.MDL_IMAGE,"%06d"%i+"@"+basisStack)
        
        # Set the class average for this class
        MDclassAveragesAux.addObject()
        MDclassAveragesAux.setValue(xmipp.MDL_IMAGE,"%06d"%0+"@"+basisStack)
        MDclassAverages.addObject()
        MDclassAverages.setValue(xmipp.MDL_IMAGE,"%06d"%idx+"@"+fnRoot+"_core.stk")
        MDaux=xmipp.MetaData()
        MDaux.readBlock(classFile,block+"_filtered")
        MDclassAverages.setValue(xmipp.MDL_IMAGE_CLASS_COUNT,int(MDaux.size()))
        classIdx[block]=idx
        idx+=1
    MDbasis.write(fnRoot+"_core_basis.sel")
    MDclassAverages.write(fnRoot+"_core.sel")
    MDclassAveragesAux.write(fnRoot+"_coreAux.sel")
    os.system("xmipp_convert_image -i "+fnRoot+"_coreAux.sel -o "+fnRoot+"_core.stk -v 0; "+\
        "rm -f "+fnRoot+"_coreAux.sel")
    os.system("xmipp_convert_image -i "+fnRoot+"_core_basis.sel -o "+fnRoot+"_core_basis.stk -v 0; "+\
        "rm -f "+fnRoot+"_core_basis.sel")
    classDir,classRoot=os.path.split(fnRoot)
    fnDirCoreClasses=classDir+"/class_core_classes"
    if os.path.exists(fnDirCoreClasses):
        os.system('rm -rf '+fnDirCoreClasses)
    os.mkdir(fnDirCoreClasses)

    for id in coreInformation:
        block=id[0]
        alignedStack=id[1]
        PCAoutput=id[3]

        # Read the PCA aligned images and build a dictionary
        MDpca=xmipp.MetaData(PCAoutput)
        correspondances={}
        for aux in MDpca:
            correspondances[MDpca.getValue(xmipp.MDL_IMAGE_ORIGINAL)]=\
                MDpca.getValue(xmipp.MDL_IMAGE)
        
        # Make the correspondance between images in the alignedStack and in the core
        MDcore=xmipp.MetaData()
        MDcore.readBlock(classFile,block+"_filtered")
        MDcoreFull=xmipp.MetaData()
        MDcoreAux=xmipp.MetaData()
        idx=0
        fnClassStack=classDir+"/class_core_classes/class_core_%05d"%classIdx[block]+".stk"
        for aux in MDcore:
            fnImageStackAligned=MDcore.getValue(xmipp.MDL_IMAGE)
            fnImageClassStack="%06d"%idx+"@"+fnClassStack
            fnImagePCAStack=correspondances[fnImageStackAligned]
            idx+=1
            
            MDcoreAux.addObject()
            MDcoreAux.setValue(xmipp.MDL_IMAGE,fnImagePCAStack)

            MDcoreFull.addObject()
            MDcoreFull.setValue(xmipp.MDL_IMAGE,fnImageClassStack)
            MDcoreFull.setValue(xmipp.MDL_ZSCORE,MDcore.getValue(xmipp.MDL_ZSCORE))
            MDcoreFull.setValue(xmipp.MDL_IMAGE_ORIGINAL,fnImageStackAligned)
            MDcoreFull.setValue(xmipp.MDL_IMAGE_CLASS,"%06d"%classIdx[block]+"@"+fnRoot+"_core.stk")
            
        MDcoreAux.write(fnRoot+"_inter.sel")
        if MDcoreAux.size()>0:
            os.system("xmipp_convert_image -i "+fnRoot+"_inter.sel -o "+fnClassStack+" -v 0")
        os.system("rm -f "+fnRoot+"_inter.sel "+PCAoutput)
        MDcoreFull.writeBlock(fnRoot+"_core.sel",block+"_filtered")
    os.system("rm -rf "+classDir+"/Basis")

# --------------------------------------------------------------------
def summary(fnRoot,Nlevels):
    print "Summary ---------------------------------------------------"
    fh=open(fnRoot+"_analysis_summary.txt",'w')
    for level in range(1,Nlevels):
        fh.write("Level: "+str(level)+'\n')
        classFile=fnRoot+"_level_%02d.sel"%level
        stackFile=classFile.replace(".sel",".stk")
        Nclasses=xmipp.SingleImgSize(stackFile)[3]
        totalSum=0
        finalSum=0
        goodClasses=0;
        for nclass in range(0,Nclasses):
            blockName="class_%06d"%nclass
            MD=xmipp.MetaData()
            MD.readBlock(classFile,blockName)
            Ninitial=MD.size()
            MD.readBlock(classFile,blockName+"_filtered")
            Nintermediate=MD.size()
            MD.readBlock(classFile,blockName+"_core")
            Nfinal=MD.size()
            if Nfinal>0:
                goodClasses=goodClasses+1
            fh.write("     Class: "+str(nclass)+": "+str(Ninitial)+\
                " -> "+str(Nintermediate)+' (filtered) -> '+str(Nfinal)+' (core)\n')
            totalSum=totalSum+Ninitial
            finalSum=finalSum+Nfinal
        fh.write("     There are "+str(goodClasses)+" good classes at this level\n")
        fh.write( "     The number of images represented has shrinked from "+\
            str(totalSum)+" to "+str(finalSum)+'\n')
    
    # After filtering the cores
    fh.write("Final classification --------------------------------\n")
    classFile=fnRoot+".sel"
    stackFile=classFile.replace(".sel",".stk")
    Nclasses=xmipp.SingleImgSize(stackFile)[3]
    totalSum=0
    finalSum=0
    goodClasses=0;
    for nclass in range(0,Nclasses):
        blockName="class_%06d"%nclass
        MD=xmipp.MetaData()
        MD.readBlock(classFile,blockName)
        Ninitial=MD.size()
        MD.readBlock(classFile,blockName+"_filtered")
        Nintermediate1=MD.size()
        MD.readBlock(classFile,blockName+"_core")
        Nintermediate2=MD.size()
        MD.readBlock(classFile,blockName+"_core_filtered")
        Nfinal=MD.size()
        if Nfinal>0:
            goodClasses=goodClasses+1
        fh.write("     Class: "+str(nclass)+": "+str(Ninitial)+\
            " -> "+str(Nintermediate1)+' (filtered) '+
            " -> "+str(Nintermediate2)+' (core) '+
            ' -> '+str(Nfinal)+' (core filtered)\n')
        totalSum=totalSum+Ninitial
        finalSum=finalSum+Nfinal
    fh.write("     There are finally "+str(goodClasses)+" good classes\n")
    fh.write( "     The number of images represented has shrinked from "+\
        str(totalSum)+" to "+str(finalSum)+'\n')
    
    fh.close()

# --------------------------------------------------------------------
def sortImages(fnRoot,Nproc,systemFlavour):
    print "Sorting images --------------------------------------------"
    import launch_job
    os.system(launch_job.buildCommand('xmipp_sort_images',
               '-i '+fnRoot+"_core.sel --oroot "+fnRoot+"_core_sorted",
               True,Nproc,1,systemFlavour))

# --------------------------------------------------------------------

def oldMain():
    if not sys.argv[1:] or len(sys.argv)<=5:
        print "Usage: xmipp_classify_CL2D_core_analysis <rootname> <thGoodClass> <thJunkZscore> <thPCAZscore> <Nproc> <system-flavour=''>"
        sys.exit()
    args = sys.argv[1:]
    fnRoot=args[0]
    thGoodClass=float(args[1])
    thGoodClass=thGoodClass/100
    thZscore=float(args[2])
    thPCAZscore=float(args[3])
    Nproc=int(args[4])
    if Nproc>1:
        Nproc=Nproc+1
    if len(args)>=6:
        systemFlavour=args[5]
    else:
        systemFlavour=''

    # Analyze the clusters as they currently are
    classFiles=sorted(glob.glob(fnRoot+"_level_??.sel"));
    classBlockRefs=[]
    for classFile in classFiles:
        stackFile=classFile.replace(".sel",".stk")
        numberOfCodes=xmipp.SingleImgSize(stackFile)[3]
        for codeNumber in range(numberOfCodes):
            blockName="class_%06d"%codeNumber
            reference=str(codeNumber)+"@"+stackFile
            classBlockRefs.append((blockName,classFile,reference))
    removeOutliersWithinCluster(Nproc,fnRoot,classBlockRefs,thZscore,thPCAZscore,
        False,systemFlavour);

    # Compute the core of each level
    Nlevels=len(classFiles);
    computeCores(1,fnRoot,Nlevels);
    
    # Purify cores
    classFile=fnRoot+".sel";
    classBlockRefs=[]
    stackFile=classFile.replace(".sel",".stk")
    numberOfCodes=xmipp.SingleImgSize(stackFile)[3]
    for codeNumber in range(numberOfCodes):
        blockName="class_%06d"%codeNumber+"_core"
        reference=str(codeNumber)+"@"+stackFile
        classBlockRefs.append((blockName,classFile,reference))
    coreInformation=removeOutliersWithinCluster(Nproc,fnRoot,classBlockRefs,
        thZscore,thPCAZscore,True,systemFlavour);

    # Gather core information
    gatherCores(fnRoot,classFile,coreInformation)

    # Summary
    summary(fnRoot,Nlevels)
    
    # Sort images
    sortImages(fnRoot,Nproc,systemFlavour)
   

from protlib_xmipp import XmippScript

class ScriptCL2DCoreAnalysis(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)

    def defineParams(self):
        self.addUsageLine("This programs needs to be re-written with new libs utilities")
        
    def run(self):
        print "NOT-IMPLEMENTED: This programs needs to be re-written with new libs utilities"
         
if __name__ == "__main__":
    ScriptCL2DCoreAnalysis().run()
