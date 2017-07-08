// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooDSCBDict

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
#include "RooDSCB.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RooDSCB(void *p = 0);
   static void *newArray_RooDSCB(Long_t size, void *p);
   static void delete_RooDSCB(void *p);
   static void deleteArray_RooDSCB(void *p);
   static void destruct_RooDSCB(void *p);
   static void streamer_RooDSCB(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooDSCB*)
   {
      ::RooDSCB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooDSCB >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooDSCB", ::RooDSCB::Class_Version(), "RooDSCB.h", 8,
                  typeid(::RooDSCB), DefineBehavior(ptr, ptr),
                  &::RooDSCB::Dictionary, isa_proxy, 16,
                  sizeof(::RooDSCB) );
      instance.SetNew(&new_RooDSCB);
      instance.SetNewArray(&newArray_RooDSCB);
      instance.SetDelete(&delete_RooDSCB);
      instance.SetDeleteArray(&deleteArray_RooDSCB);
      instance.SetDestructor(&destruct_RooDSCB);
      instance.SetStreamerFunc(&streamer_RooDSCB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooDSCB*)
   {
      return GenerateInitInstanceLocal((::RooDSCB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooDSCB*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooDSCB::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooDSCB::Class_Name()
{
   return "RooDSCB";
}

//______________________________________________________________________________
const char *RooDSCB::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooDSCB*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooDSCB::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooDSCB*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooDSCB::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooDSCB*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooDSCB::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooDSCB*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooDSCB::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooDSCB.

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
      R__b.CheckByteCount(R__s, R__c, RooDSCB::IsA());
   } else {
      R__c = R__b.WriteVersion(RooDSCB::IsA(), kTRUE);
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
   static void *new_RooDSCB(void *p) {
      return  p ? new(p) ::RooDSCB : new ::RooDSCB;
   }
   static void *newArray_RooDSCB(Long_t nElements, void *p) {
      return p ? new(p) ::RooDSCB[nElements] : new ::RooDSCB[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooDSCB(void *p) {
      delete ((::RooDSCB*)p);
   }
   static void deleteArray_RooDSCB(void *p) {
      delete [] ((::RooDSCB*)p);
   }
   static void destruct_RooDSCB(void *p) {
      typedef ::RooDSCB current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooDSCB(TBuffer &buf, void *obj) {
      ((::RooDSCB*)obj)->::RooDSCB::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooDSCB

namespace {
  void TriggerDictionaryInitialization_RooDSCBDict_Impl() {
    static const char* headers[] = {
"RooDSCB.h",
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
class __attribute__((annotate("$clingAutoload$RooDSCB.h")))  RooDSCB;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "RooDSCB.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooDSCB", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooDSCBDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooDSCBDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooDSCBDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooDSCBDict() {
  TriggerDictionaryInitialization_RooDSCBDict_Impl();
}
