#include "R3BCalifaContFact.h"

#include "FairLogger.h"
#include "FairParAsciiFileIo.h"
#include "FairParRootFileIo.h"
#include "FairRuntimeDb.h"
#include "R3BCalifaCrystalCalPar.h"

#include "TClass.h"

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

ClassImp(R3BCalifaContFact)

    static R3BCalifaContFact gR3BCalifaContFact;

R3BCalifaContFact::R3BCalifaContFact()
{
    /** Constructor (called when the library is loaded) */
    fName = "R3BCalifaContFact";
    fTitle = "Factory for parameter containers in libR3BCalifa";
    setAllContainers();
    FairRuntimeDb::instance()->addContFactory(this);
}

void R3BCalifaContFact::setAllContainers()
{
    /** Creates the Container objects with all accepted contexts and adds them
        to the list of containers for the CALIFA library.
    */

    FairContainer* p1 =
        new FairContainer("califaCrystalCalPar", "Califa Calibration Parameters", "CalifaCalParContext");
    p1->addContext("CalifaCalParContext");

    containers->Add(p1);
}

FairParSet* R3BCalifaContFact::createContainer(FairContainer* c)
{
    /** Calls the constructor of the corresponding parameter container.
     * For an actual context, which is not an empty string and not the default context
     * of this container, the name is concatinated with the context.
     */

    const char* name = c->GetName();
    LOG(INFO) << "R3BCalifaContFact: Create container name: " << name;
    FairParSet* p = 0;
    if (strcmp(name, "califaCrystalCalPar") == 0)
    {
        p = new R3BCalifaCrystalCalPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
    }
    return p;
}

void R3BCalifaContFact::activateParIo(FairParIo* io)
{
    // activates the input/output class for the parameters
    // needed by the Sts
    /*
      if (strcmp(io->IsA()->GetName(),"FairParRootFileIo")==0) {
      R3BCalifaParRootFileIo* p=new R3BCalifaParRootFileIo(((FairParRootFileIo*)io)->getParRootFile());
      io->setDetParIo(p);
      }
      if (strcmp(io->IsA()->GetName(),"FairParAsciiFileIo")==0) {
      R3BCalifaParAsciiFileIo* p=new R3BCalifaParAsciiFileIo(((FairParAsciiFileIo*)io)->getFile());
      io->setDetParIo(p);
      }
    */
}
