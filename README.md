### EIC Jet Skeleton ###

##### About #####

Demonstrate usage of smeared trees from eic-smear in some more detail, with a simple jet analysis.
The installation is intended to be independent of cmake etc., with a simple Makefile.


Author / Contact:
* Kolja Kauder <kkauder@bnl.gov>

##### Prerequisites #####

* An existing eic-smear installation, assumed to be pointed to by $EICDIRECTORY.
* A ROOT installation, assumed to be pointed to by $ROOTSYS
* Truth and smeared MC trees, in this documentation assumed to be named "truth.root" and "truth.smeared.root".
An easy way to obtain some with a standard installation is:
```
cp $EICDIRECTORY/PACKAGES/eic-smear/tests/ep_hiQ2.20x250.small.txt truth.txt
$EICDIRECTORY/PACKAGES/eic-smear/build/tests/qaplots -i truth.txt
```
* A fastjet3 installation, pointed to by $FASTJETDIR
* Optionally, an fjcontrib installation. The Makefile will check for that.

##### Build #####
```
make
```



