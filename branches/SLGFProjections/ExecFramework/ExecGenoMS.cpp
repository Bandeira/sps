//Header Include
#include "ExecGenoMS.h"

// Module Includes
#include "ExecBase.h"
#include "utils.h"
#include "Logger.h"

// System Includes
#include <stdio.h>
#include <string.h>


using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ExecGenoMS::ExecGenoMS(void)
{
  m_name = "ExecGenoMS";
  m_type = "ExecGenoMS";
  m_dataOwned = true;
  m_isValid = false;
}

// -------------------------------------------------------------------------
ExecGenoMS::ExecGenoMS(const ParameterList & params) 
{
  m_name = "ExecGenoMS";
  m_type = "ExecGenoMS";
  m_dataOwned = true;
  m_isValid = false;
  m_params = params;
}

// -------------------------------------------------------------------------
ExecGenoMS::~ExecGenoMS(void)
{
  if(m_dataOwned)
    {
      //delete shit
    }

}

// -------------------------------------------------------------------------
bool ExecGenoMS::isValid(void) const
{
  return m_isValid;
}

// -------------------------------------------------------------------------
string ExecGenoMS::getName(void) const
{
  return m_name;
}

// -------------------------------------------------------------------------
string ExecGenoMS::getType(void) const
{
  return m_type;
}

// -------------------------------------------------------------------------
void ExecGenoMS::setName(std::string name)
{
  m_name = name;
}

// -------------------------------------------------------------------------

//Full validation of the parameters happens within GenoMS code
bool ExecGenoMS::validateParams(std::string & error)
{

  m_isValid = false;

  VALIDATE_PARAM_EXIST("OUTPUT_FILE");
  VALIDATE_PARAM_EXIST("EXE_DIR");

  m_isValid = true;

 
  return true;
}



ExecBase * ExecGenoMS::clone(const ParameterList & input_params) const
{
  return new ExecGenoMS(input_params);
}

bool ExecGenoMS::loadInputData(void)
{
  return true;
}

 
bool ExecGenoMS::saveOutputData(void)
{
  return true;
}

bool ExecGenoMS::saveInputData(std::vector<std::string> & filenames)
{
  return true;
}

bool ExecGenoMS::loadOutputData(void)
{
  return true;
}

vector<ExecBase *> const & ExecGenoMS::split(int numSplit)
{
  m_subModules.resize(0);
  return m_subModules;
}

bool ExecGenoMS::merge(void)
{
}



bool ExecGenoMS::invoke(void)
{

  DEBUG_MSG("Entering ExecGenoMS::invoke()");
  
  std::string exeDir = m_params.getValue("EXE_DIR");
  rtrim(exeDir);

  int ret = callGenoMS(exeDir);
  if(ret != 0)
    {
      ERROR_MSG("Error executing GenoMS!");
      return false;
    }


  return true;

}/*
JNIEnv* create_vm(JavaVM ** jvm, string & exeDir) {
    
  JNIEnv *env;
  int ret;
  string classPathStr = exeDir + PATH_SEPARATOR + "GenoMS.jar";


#ifdef JNI_VERSION_1_2

  JavaVMInitArgs vm_args;
  JavaVMOption options[4]; 
  vm_args.version = 0x00010002;
  
  //Path to the java source code     
  options[0].optionString = strdup(("-Djava.class.path=/usr/lib/jvm/jre-1.6.0-openjdk.x86_64/lib/rt.jar:/usr/lib/jvm/jre-1.6.0-openjdk.x86_64/lib:" + classPathStr + ":/usr/java/default/jre/lib").c_str());
  
  options[1].optionString = "-verbose:gc,jni,class";
  options[2].optionString = "-Xmx2000M";
  options[3].optionString = "-Djava.home=/usr/java/default/jre";
  
  vm_args.nOptions = 4;
  vm_args.options = options;
  vm_args.ignoreUnrecognized = JNI_FALSE;

  printf("args are set\n");
  ret = JNI_CreateJavaVM(jvm, (void**)&env, &vm_args);
  printf("created jvm\n");
#else 
  JDK1_1InitArgs vm_args;
  char classpath[1024];

  vm_args.version = 0x00010001;

  sprintf(classpath,"%s%c%s",vm_args.classpath, ":", classPathStr.c_str());
  vm_args.classpath = classpath;

  ret = JNI_CreateJavaVM(jvm, &env,*vm_args);

#endif 

  if(ret < 0)
    printf("\nUnable to Launch JVM: %d\n",ret);       
  else
    printf("Successful JVM Launch\n");
  return env;
}

jsize ExecGenoMS::getNumArgs()
{
 
  jsize ret = 8;
  
  //GenoMSCommand += " -i genoMS.params";
  
  //GenoMSCommand += " -o ";
  //GenoMSCommand += m_params.getValue("OUTPUT_FILE");

  //GenoMSCommand += " -r ";
  //GenoMSCommand += exeDir;
  //GenoMSCommand += "/DBs_GenoMS";
  
  //GenoMSCommand += " -e ";
  //GenoMSCommand += exeDir;


  if(m_params.exists("PROJECT_DIR"))
    { 
      ret += 2;
      //GenoMSCommand += " -k ";
      //GenoMSCommand += m_params.getValue("PROJECT_DIR");
    }

  if(m_params.exists("RUN_DBSEARCH"))
    ret += 1;
    //GenoMSCommand += " -x";
  
  if(m_params.exists("PEAK_PENALTY"))
    ret += 1; //GenoMSCommand += " -p";
 
  if(m_params.exists("ALL2ALL_SIMILARITY"))
    ret += 1; //GenoMSCommand += " -s";
  
  if(m_params.exists("HMM_LATE_ADD"))
    ret += 1; //GenoMSCommand += " -a";
  
  if(m_params.exists("FDR_CUTOFF"))
    {
      ret += 2;
      //GenoMSCommand += " -w ";
      //GenoMSCommand += m_params.getValue("FDR_CUTOFF");
    }

  if(m_params.exists("MUTATION_MODE"))
    ret += 1;//GenoMSCommand += " -f";

  if(m_params.exists("GENERATE_REPORTS"))
    {
      ret += 1;//GenoMSCommand += " -q";
    }
  if(m_params.exists("LOG_FILE"))
    {
      ret += 2;
      //GenoMSCommand += " -l ";
      //GenoMSCommand += m_params.getValue("LOG_FILE");
    }
  return ret;
}
*/
/*
int ExecGenoMS::callGenoMS(string & exeDir)
{
  JNIEnv *env;
  JavaVM * jvm;
  env = create_vm(&jvm, exeDir);
  if (env == NULL)
    return 25;  
  

  jclass cls = env->FindClass("AntibodyDriver"); 
  jmethodID mid = env->GetStaticMethodID(cls, "main", "([Ljava/lang/String;)V");

  jclass str = env->FindClass("java/lang/String");
  jsize num_args = getNumArgs();
  jobjectArray jargs = env->NewObjectArray(num_args, str, NULL);

  int currArg = 0;
  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-i"));
  currArg += 1;
  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("genoMS.params"));
  currArg += 1;

  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-o"));
  currArg += 1;
  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF(m_params.getValue("OUTPUT_FILE").c_str()));
  currArg += 1;

  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-r"));
  currArg += 1;
  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF((exeDir + PATH_SEPARATOR + "DBs_GenoMS").c_str()));
  currArg += 1;

  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-e"));
  currArg += 1;
  env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF(exeDir.c_str()));
  currArg += 1;


  if(m_params.exists("PROJECT_DIR"))
    { 
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-k"));
      currArg += 1;
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF(m_params.getValue("PROJECT_DIR").c_str()));
      currArg += 1;
    }

  if(m_params.exists("RUN_DBSEARCH"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-x"));
      currArg += 1;
    }
  //GenoMSCommand += " -x";
  
  if(m_params.exists("PEAK_PENALTY"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-p"));
      currArg += 1;
    }

   if(m_params.exists("ALL2ALL_SIMILARITY"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-s"));
      currArg += 1;
    }
  //GenoMSCommand += " -s";
  
  if(m_params.exists("HMM_LATE_ADD"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-a"));
      currArg += 1;
    }
  //    GenoMSCommand += " -a";
  
  if(m_params.exists("FDR_CUTOFF"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-w"));
      currArg += 1;
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF(m_params.getValue("FDR_CUTOFF").c_str()));
      currArg += 1;
      //  GenoMSCommand += " -w ";
      //GenoMSCommand += m_params.getValue("FDR_CUTOFF");
    }

  if(m_params.exists("MUTATION_MODE"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-f"));
      currArg += 1;
    }
  //GenoMSCommand += " -f";

  if(m_params.exists("GENERATE_REPORTS"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-q"));
      currArg += 1;
      //GenoMSCommand += " -q";
    }
  if(m_params.exists("LOG_FILE"))
    {
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF("-l"));
      currArg += 1;
      env->SetObjectArrayElement(jargs, currArg, env->NewStringUTF(m_params.getValue("LOG_FILE").c_str()));
      currArg += 1;
      //GenoMSCommand += " -l ";
      //GenoMSCommand += m_params.getValue("LOG_FILE");
    }

  env->CallStaticVoidMethod(cls, mid, jargs); 

  jvm->DestroyJavaVM();
}
*/

int ExecGenoMS::callGenoMS(string & exeDir)
{
  //GenoMS command
  std::string GenoMSCommand;

  //Because I don't know how to get the name of the original config file,
  //We have to write a new one for GenoMS
  m_params.setValue("MSGFDBPATH",exeDir);
  m_params.writeToFile("genoMS.params");

  GenoMSCommand = "unset LD_LIBRARY_PATH && java -jar ";
  GenoMSCommand += exeDir;
  GenoMSCommand += "/GenoMS.jar";
  
  GenoMSCommand += " -i genoMS.params";
  
  GenoMSCommand += " -o ";
  GenoMSCommand += m_params.getValue("OUTPUT_FILE");

  GenoMSCommand += " -r ";
  GenoMSCommand += exeDir;
  GenoMSCommand += "/DBs_GenoMS";

  if(m_params.exists("PROJECT_DIR"))
    { 
      GenoMSCommand += " -k ";
      GenoMSCommand += m_params.getValue("PROJECT_DIR");
    }

  GenoMSCommand += " -e ";
  GenoMSCommand += exeDir;
  

  if(m_params.exists("RUN_DBSEARCH"))
    GenoMSCommand += " -x";
  
  if(m_params.exists("PEAK_PENALTY"))
    GenoMSCommand += " -p";
 
  if(m_params.exists("ALL2ALL_SIMILARITY"))
    GenoMSCommand += " -s";
  
  if(m_params.exists("HMM_LATE_ADD"))
    GenoMSCommand += " -a";
  
  if(m_params.exists("FDR_CUTOFF"))
    {
      GenoMSCommand += " -w ";
      GenoMSCommand += m_params.getValue("FDR_CUTOFF");
    }

  if(m_params.exists("MUTATION_MODE"))
    GenoMSCommand += " -f";

  if(m_params.exists("GENERATE_REPORTS"))
    {
      GenoMSCommand += " -q";
    }
  if(m_params.exists("LOG_FILE"))
    {
      GenoMSCommand += " -l ";
      GenoMSCommand += m_params.getValue("LOG_FILE");
    }
  
  DEBUG_MSG("call genoMS: " << GenoMSCommand);
  
  int ret = system(GenoMSCommand.c_str()); 
  DEBUG_MSG("GenoMS return value: " << ret);
  return ret;
  
}
