#include "GNNTrakerCom.h"
#include "TrakerImp.h"
#include "stdafx.h"

extern "C" GNN_Traker * createGNNTracker()
{
	return new TrakerImp;
}
