// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once


#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#define _USE_MATH_DEFINES //IN order to use Math COnstants
#include <stdio.h>
#include <tchar.h>
#include <stdlib.h> // for srand ( ) and rand ( ) 
#include <time.h> // for time ( ) and time_t 
#include <iostream>
#include <fstream> //File Streams
#include <direct.h> // for getcwd 
//#include <math.h>     // for exp(), log(), and log10()

// ###GSL Note: For the library to work in MSVC, I had to change to the Multithreaded version WinGsl_md.lib
// Also under Properties->C/C++->Code GEneration->Run Time Library Change to Multithreaded Debug
#include <WinGsl.h >

///Program Parameters

#define MAX_AFFERENTS 1202
#define MAX_SPIKES	2500 //Number of simultanuous spike that can be injecting to the neuron 
#define G_MAX		0.0151  //0.015 //The Song Conductance Max VAlue used in IFNeuron And SynapseEnsemble
#define G_INH		0.05 //0.05 //The Song INH Conductance Fixed VAlue used in IFNeuron And SynapseEnsemble
#define IFFIRERATE_PERIOD 1 ///The time over which the IFNeuron calculates the Fire rate
//#define DEBUG_LOG
#define USE_SONG_CONDUCTANCE //When Defined the simpler Song method is used to calculate gex
//#define USE_SONG_LEARNING	// Synaptic modification implemented as the double exponential rule and switch rule is ignored.
//#define VERBOSE

// TODO: reference additional headers your program requires here
