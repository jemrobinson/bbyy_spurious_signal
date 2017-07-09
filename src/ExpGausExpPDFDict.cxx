// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIExpGausExpPDFDict

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
#include "ExpGausExpPDF.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_ExpGausExpPDF(void *p = 0);
   static void *newArray_ExpGausExpPDF(Long_t size, void *p);
   static void delete_ExpGausExpPDF(void *p);
   static void deleteArray_ExpGausExpPDF(void *p);
   static void destruct_ExpGausExpPDF(void *p);
   static void streamer_ExpGausExpPDF(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ExpGausExpPDF*)
   {
      ::ExpGausExpPDF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ExpGausExpPDF >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ExpGausExpPDF", ::ExpGausExpPDF::Class_Version(), "ExpGausExpPDF.h", 8,
                  typeid(::ExpGausExpPDF), DefineBehavior(ptr, ptr),
                  &::ExpGausExpPDF::Dictionary, isa_proxy, 16,
                  sizeof(::ExpGausExpPDF) );
      instance.SetNew(&new_ExpGausExpPDF);
      instance.SetNewArray(&newArray_ExpGausExpPDF);
      instance.SetDelete(&delete_ExpGausExpPDF);
      instance.SetDeleteArray(&deleteArray_ExpGausExpPDF);
      instance.SetDestructor(&destruct_ExpGausExpPDF);
      instance.SetStreamerFunc(&streamer_ExpGausExpPDF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ExpGausExpPDF*)
   {
      return GenerateInitInstanceLocal((::ExpGausExpPDF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ExpGausExpPDF*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ExpGausExpPDF::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ExpGausExpPDF::Class_Name()
{
   return "ExpGausExpPDF";
}

//______________________________________________________________________________
const char *ExpGausExpPDF::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExpGausExpPDF*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ExpGausExpPDF::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ExpGausExpPDF*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ExpGausExpPDF::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExpGausExpPDF*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ExpGausExpPDF::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ExpGausExpPDF*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ExpGausExpPDF::Streamer(TBuffer &R__b)
{
   // Stream an object of class ExpGausExpPDF.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      m0.Streamer(R__b);
      sigma.Streamer(R__b);
      kLo.Streamer(R__b);
      kHi.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, ExpGausExpPDF::IsA());
   } else {
      R__c = R__b.WriteVersion(ExpGausExpPDF::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      m0.Streamer(R__b);
      sigma.Streamer(R__b);
      kLo.Streamer(R__b);
      kHi.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ExpGausExpPDF(void *p) {
      return  p ? new(p) ::ExpGausExpPDF : new ::ExpGausExpPDF;
   }
   static void *newArray_ExpGausExpPDF(Long_t nElements, void *p) {
      return p ? new(p) ::ExpGausExpPDF[nElements] : new ::ExpGausExpPDF[nElements];
   }
   // Wrapper around operator delete
   static void delete_ExpGausExpPDF(void *p) {
      delete ((::ExpGausExpPDF*)p);
   }
   static void deleteArray_ExpGausExpPDF(void *p) {
      delete [] ((::ExpGausExpPDF*)p);
   }
   static void destruct_ExpGausExpPDF(void *p) {
      typedef ::ExpGausExpPDF current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ExpGausExpPDF(TBuffer &buf, void *obj) {
      ((::ExpGausExpPDF*)obj)->::ExpGausExpPDF::Streamer(buf);
   }
} // end of namespace ROOT for class ::ExpGausExpPDF

namespace {
  void TriggerDictionaryInitialization_ExpGausExpPDFDict_Impl() {
    static const char* headers[] = {
"ExpGausExpPDF.h",
0
    };
    static const char* includePaths[] = {
"./include",
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
class __attribute__((annotate("$clingAutoload$ExpGausExpPDF.h")))  ExpGausExpPDF;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "ExpGausExpPDF.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ExpGausExpPDF", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ExpGausExpPDFDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ExpGausExpPDFDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ExpGausExpPDFDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ExpGausExpPDFDict() {
  TriggerDictionaryInitialization_ExpGausExpPDFDict_Impl();
}
