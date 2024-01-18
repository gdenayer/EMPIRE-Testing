//EMPIRE_API.c
#include "empire.h"
#include "EMPIRE_API.h"

EMPIRE *empire;

void EMPIRE_Init(){
  empire = new EMPIRE();
  empire->Init();
}

void EMPIRE_Send(){
      empire->Send();
}
void EMPIRE_Receive(){
      empire->Receive();
}
void EMPIRE_Disconnect(){
  empire->Disconnect();
  delete empire;
}