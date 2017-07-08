// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RooExpGausExpDict

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
#include "RooExpGausExp.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RooExpGausExp(void *p = 0);
   static void *newArray_RooExpGausExp(Long_t size, void *p);
   static void delete_RooExpGausExp(void *p);
   static void deleteArray_RooExpGausExp(void *p);
   static void destruct_RooExpGausExp(void *p);
   static void streamer_RooExpGausExp(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooExpGausExp*)
   {
      ::RooExpGausExp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooExpGausExp >(0);
      static ::ROOT::TGenericClassInfo
         instance("RooExpGausExp", ::RooExpGausExp::Class_Version(), "RooExpGausExp.h", 8,
                  typeid(::RooExpGausExp), DefineBehavior(ptr, ptr),
                  &::RooExpGausExp::Dictionary, isa_proxy, 16,
                  sizeof(::RooExpGausExp) );
      instance.SetNew(&new_RooExpGausExp);
      instance.SetNewArray(&newArray_RooExpGausExp);
      instance.SetDelete(&delete_RooExpGausExp);
      instance.SetDeleteArray(&deleteArray_RooExpGausExp);
      instance.SetDestructor(&destruct_RooExpGausExp);
      instance.SetStreamerFunc(&streamer_RooExpGausExp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooExpGausExp*)
   {
      return GenerateInitInstanceLocal((::RooExpGausExp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooExpGausExp*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooExpGausExp::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooExpGausExp::Class_Name()
{
   return "RooExpGausExp";
}

//______________________________________________________________________________
const char *RooExpGausExp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooExpGausExp*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooExpGausExp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooExpGausExp*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooExpGausExp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooExpGausExp*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooExpGausExp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooExpGausExp*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooExpGausExp::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooExpGausExp.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      m0.Streamer(R__b);
      sigma.Streamer(R__b);
      kLo.Streamer(R__b);
      kHi.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooExpGausExp::IsA());
   } else {
      R__c = R__b.WriteVersion(RooExpGausExp::IsA(), kTRUE);
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
   static void *new_RooExpGausExp(void *p) {
      return  p ? new(p) ::RooExpGausExp : new ::RooExpGausExp;
   }
   static void *newArray_RooExpGausExp(Long_t nElements, void *p) {
      return p ? new(p) ::RooExpGausExp[nElements] : new ::RooExpGausExp[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooExpGausExp(void *p) {
      delete ((::RooExpGausExp*)p);
   }
   static void deleteArray_RooExpGausExp(void *p) {
      delete [] ((::RooExpGausExp*)p);
   }
   static void destruct_RooExpGausExp(void *p) {
      typedef ::RooExpGausExp current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooExpGausExp(TBuffer &buf, void *obj) {
      ((::RooExpGausExp*)obj)->::RooExpGausExp::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooExpGausExp

namespace {
  void TriggerDictionaryInitialization_RooExpGausExpDict_Impl() {
    static const char* headers[] = {
"RooExpGausExp.h",
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
class __attribute__((annotate("$clingAutoload$RooExpGausExp.h")))  RooExpGausExp;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "RooExpGausExp.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooExpGausExp", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RooExpGausExpDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RooExpGausExpDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RooExpGausExpDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RooExpGausExpDict() {
  TriggerDictionaryInitialization_RooExpGausExpDict_Impl();
}
