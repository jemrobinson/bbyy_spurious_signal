// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DoubleSidedCrystalBallPDFDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "DoubleSidedCrystalBallPDF.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_DoubleSidedCrystalBallPDF(void *p = 0);
   static void *newArray_DoubleSidedCrystalBallPDF(Long_t size, void *p);
   static void delete_DoubleSidedCrystalBallPDF(void *p);
   static void deleteArray_DoubleSidedCrystalBallPDF(void *p);
   static void destruct_DoubleSidedCrystalBallPDF(void *p);
   static void streamer_DoubleSidedCrystalBallPDF(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DoubleSidedCrystalBallPDF*)
   {
      ::DoubleSidedCrystalBallPDF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DoubleSidedCrystalBallPDF >(0);
      static ::ROOT::TGenericClassInfo
         instance("DoubleSidedCrystalBallPDF", ::DoubleSidedCrystalBallPDF::Class_Version(), "DoubleSidedCrystalBallPDF.h", 8,
                  typeid(::DoubleSidedCrystalBallPDF), DefineBehavior(ptr, ptr),
                  &::DoubleSidedCrystalBallPDF::Dictionary, isa_proxy, 16,
                  sizeof(::DoubleSidedCrystalBallPDF) );
      instance.SetNew(&new_DoubleSidedCrystalBallPDF);
      instance.SetNewArray(&newArray_DoubleSidedCrystalBallPDF);
      instance.SetDelete(&delete_DoubleSidedCrystalBallPDF);
      instance.SetDeleteArray(&deleteArray_DoubleSidedCrystalBallPDF);
      instance.SetDestructor(&destruct_DoubleSidedCrystalBallPDF);
      instance.SetStreamerFunc(&streamer_DoubleSidedCrystalBallPDF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DoubleSidedCrystalBallPDF*)
   {
      return GenerateInitInstanceLocal((::DoubleSidedCrystalBallPDF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DoubleSidedCrystalBallPDF*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr DoubleSidedCrystalBallPDF::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DoubleSidedCrystalBallPDF::Class_Name()
{
   return "DoubleSidedCrystalBallPDF";
}

//______________________________________________________________________________
const char *DoubleSidedCrystalBallPDF::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DoubleSidedCrystalBallPDF*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DoubleSidedCrystalBallPDF::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DoubleSidedCrystalBallPDF*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DoubleSidedCrystalBallPDF::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DoubleSidedCrystalBallPDF*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DoubleSidedCrystalBallPDF::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DoubleSidedCrystalBallPDF*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void DoubleSidedCrystalBallPDF::Streamer(TBuffer &R__b)
{
   // Stream an object of class DoubleSidedCrystalBallPDF.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      m0.Streamer(R__b);
      sigma.Streamer(R__b);
      alphaLo.Streamer(R__b);
      nLo.Streamer(R__b);
      alphaHi.Streamer(R__b);
      nHi.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, DoubleSidedCrystalBallPDF::IsA());
   } else {
      R__c = R__b.WriteVersion(DoubleSidedCrystalBallPDF::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      m0.Streamer(R__b);
      sigma.Streamer(R__b);
      alphaLo.Streamer(R__b);
      nLo.Streamer(R__b);
      alphaHi.Streamer(R__b);
      nHi.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DoubleSidedCrystalBallPDF(void *p) {
      return  p ? new(p) ::DoubleSidedCrystalBallPDF : new ::DoubleSidedCrystalBallPDF;
   }
   static void *newArray_DoubleSidedCrystalBallPDF(Long_t nElements, void *p) {
      return p ? new(p) ::DoubleSidedCrystalBallPDF[nElements] : new ::DoubleSidedCrystalBallPDF[nElements];
   }
   // Wrapper around operator delete
   static void delete_DoubleSidedCrystalBallPDF(void *p) {
      delete ((::DoubleSidedCrystalBallPDF*)p);
   }
   static void deleteArray_DoubleSidedCrystalBallPDF(void *p) {
      delete [] ((::DoubleSidedCrystalBallPDF*)p);
   }
   static void destruct_DoubleSidedCrystalBallPDF(void *p) {
      typedef ::DoubleSidedCrystalBallPDF current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DoubleSidedCrystalBallPDF(TBuffer &buf, void *obj) {
      ((::DoubleSidedCrystalBallPDF*)obj)->::DoubleSidedCrystalBallPDF::Streamer(buf);
   }
} // end of namespace ROOT for class ::DoubleSidedCrystalBallPDF

namespace {
  void TriggerDictionaryInitialization_DoubleSidedCrystalBallPDFDict_Impl() {
    static const char* headers[] = {
"DoubleSidedCrystalBallPDF.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/include",
"/afs/cern.ch/user/j/jrobinso/HGamma/SpuriousSignal/",
0
    };
    static const char* fwdDeclCode =
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$DoubleSidedCrystalBallPDF.h")))  DoubleSidedCrystalBallPDF;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "DoubleSidedCrystalBallPDF.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DoubleSidedCrystalBallPDF", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DoubleSidedCrystalBallPDFDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DoubleSidedCrystalBallPDFDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DoubleSidedCrystalBallPDFDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DoubleSidedCrystalBallPDFDict() {
  TriggerDictionaryInitialization_DoubleSidedCrystalBallPDFDict_Impl();
}
