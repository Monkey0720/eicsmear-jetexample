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
If no smearing device accepted a particular particle, the corresponding entry will be a null pointer (e.g. outside acceptance, or not a final particle).
A more subtle situation arises when partial information is smeared, such as energy. In such a case, the remaining observables in the ```ParticleMCS``` object will be set to 0. Consequences:

* You should check for "true" 0 using something like ```std::abs(P)<1e-7``` to distinguish between small values and unmeasured values.
* That is especially problematic for angle values (phi, theta) and part of the reason why a detector parameterization should always include some phi, theta device for all accepted particles.
* If only partial information is available but a full 4-vector is needed, choices need to be made.

###### Possible combinations ######
* Pointer is non-zero, yet ```P=0 and E=0```.
This can happen when e. g. a low-p particle gets smeared below 0. We cannot really distinguish it from a particle that wasn't measured in the first place, but neither should we.
* Everything is measured, ```P>0 and E>0```. This is often the default case for charged particles, and the examples will use the 4-vector (E,P). Nevertheless, an analyzer could use outside knowledge and for example replace an E measurement in a region of known bad resolution with a value calculated from P and PID information.
* Tracker only, ```P>0 and E=0```. This could indicate a small E value smeared to below 0, or a region not covered by calorimetry.
Here, a real-world analyzer would know which calorimetry is missing and thus be able to narrow down the particle's genre to hadronic or electromagnetic. The framework cannot infer this information. Nevertheless, it is not unreasonable to assume that in fact HCal is missing, since an EIC detector would have large EMCal coverage. This is the approach in ```jetspectra.cxx```. In ```extendedjetexample.cxx```, we instead demonstrate how to get the genre information from the truth level. In either case, this is all that can be known without additional PID measurements, so the mass to use for an unknown particle needs to be a decision by the experimenter. Pion mass for charged hadrons and electron mass or zero for a charged electromagnetic particle are typical assumptions
* Tracker only, ```P=0 and E>0```.
The most difficult case, and best handled by acquiring the genre information from the truth level as above; a real measurement would be clearly in an HCal or an EMCal. But in either case, a calorimeter hit without a pointing track or PID information requires thoughtful handling to decide whether this is most likely a neutral hadron (neutron? K0L?) or a tracking inefficiency. This would also lead to wildly different assumptions about what the phi, theta position in the calorimeter means for the particle's true momentum (which is always smeared just marginally from the truth value). The best approach would in most cases be to reject this particle entirely, but some other options are shown in ```extendedjetexample.cxx```.

###### Other features ######
* Eic-smear does not provide tracking inefficiency. Current parameterizations also do not have pT-dependent acceptance. These depend heavily on geometry, B-field, reconstruction algorithms, quality cuts, etc, and may partially or fully be in future framework releases as more detailed concepts and standardized tracking assumptions coalesce.
In the meantime, analyzers will nevertheless want to make afterburner assumptions. An example on how to do this is provided in the form of a (completely arbitrary!) toy example:
```
eff = new TF1("eff","(x>[2]) * [0]*TMath::Erf(x-[1])",0, 100);
// mostly 99%, dropping toward small pT, sharp cutoff at 0.2
eff->SetParameters (0.99,-0.8, 0.2);
```

* Particle-level truth-smeared matching is easy via the friend mechanism, ```extendedjetexample.cxx``` has a simple version of jet-level matching.

* If fjcontrib is installed, ```extendedjetexample.cxx``` also performs a sub-jet analysis, mostly to demonstrate linking to additional external libraries.
