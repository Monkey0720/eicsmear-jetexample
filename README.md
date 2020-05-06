### EIC Jet Analysis Skeleton ###

##### About #####

Demonstrate usage of smeared trees from eic-smear in some more detail, with a simple jet analysis.
The installation is intended to be independent of cmake etc., with a short Makefile.

A jet analysis was chosen because it highlights many of the subtleties and idiosyncrasies of working with eic-smear, but the examples here are intended for general use and an appropriate starting point for most analyses.


##### Author / Contact:
* Kolja Kauder <kkauder@bnl.gov>

Please do not hesitate to contact me. If something is unclear to you, chances are it's unclear to others as well!

##### Prerequisites #####

* An existing eic-smear installation, assumed to be pointed to by $EICDIRECTORY.
* A ROOT installation, assumed to be pointed to by $ROOTSYS
* Truth and smeared MC trees, in this documentation assumed to be named "truth.root" and "truth.smeared.root".
An easy way to obtain such trees with a standard eic-smear installation is:
```
cp $EICDIRECTORY/PACKAGES/eic-smear/tests/ep_hiQ2.20x250.small.txt truth.txt
$EICDIRECTORY/PACKAGES/eic-smear/build/tests/qaplots -i truth.txt
```
* A FastJet3 installation, pointed to by $FASTJETDIR
* Optionally, an fjcontrib installation. The Makefile will check for that.

##### Build #####
```
make
```
This will generate two executables,
```
./bin/jetspectra
./bin/extendedjetexample
```

##### Discussion #####
The README file from https://github.com/eic/eic-smear provides instructions on how to generate MC data and smear it using provided example detector parameterizations, as well as ways how to modify such parameterizations.
There are also basic instructions on contents and usage of the smeared tree. The purpose of this example is to demonstrate usage in more practical detail.

It is strongly encouraged to not treat the examples as black boxes where all that's necessary is plugging in a few lines. Take note of the comments and the choices made.

###### Stand-alone or connected to truth level? ######
The smeared tree can be used by itself, quite similar to the way an experimentalist would treat real data. This is shown in ```src/jetspectra.cxx```. In many cases however, you will want to make comparisons at the particle level between truth level and its smeared counterpart. This is easily done using ROOT's friend mechanism and shown in ```src/extendedjetexample.cxx```. As discussed below, careful usage of truth level information may also be necessary to compensate for the fact that eic-smear doesn't differentiate between energy smeared in HCal or EMCal.

###### General paradigm: "zero means not measured" ######
