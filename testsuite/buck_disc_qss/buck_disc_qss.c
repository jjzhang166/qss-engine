#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __QSSIntegrator_1_p[4];
double __QSSIntegrator_1_x0 = 0;
double __QSSIntegrator_2_p[4];
double __QSSIntegrator_2_x0 = 0;
double __WSum_3_p[9];
double __WSum_3_w[2];
double __WSum_11_p[9];
double __WSum_11_w[2];
double __qss_switch_12_p[1];
double __qss_switch_12_level = 0;
double __qss_switch_13_p[1];
double __qss_switch_13_level = 0;
double __Constant_14_p[1];
double __Constant_14_k = 0;
double __Constant_15_p[1];
double __Constant_15_k = 0;
double __square_sci_16_p[3];
double __square_sci_16_amplitude = 0;
double __square_sci_16_freq = 0;
double __square_sci_16_DC = 0;
double __hysteretic_18_p[4];
double __hysteretic_18_xl = 0;
double __hysteretic_18_xu = 0;
double __hysteretic_18_yl = 0;
double __hysteretic_18_yu = 0;
double __qss_switch_19_p[1];
double __qss_switch_19_level = 0;
double __WSum_20_p[9];
double __WSum_20_w[2];
double __WSum_21_p[9];
double __WSum_21_w[2];

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = TRUE;
	 settings->method = 6;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	switch(i)
	{
		case 0:
			alg[16] = x[0];
			alg[20] = x[4];
			alg[76] = alg[16]*__WSum_20_w[0]+alg[20]*__WSum_20_w[1];
			alg[80] = alg[76];
			dx[1] = alg[80];
			return;
		case 1:
			alg[0] = __Constant_14_k;
			alg[4] = 100000.0;
			alg[24] = alg[4];
			alg[32] = alg[0];
			alg[36] = alg[0];
			alg[44] = alg[4];
			alg[48] = alg[36]*d[(0)]+alg[44]*(1.0-d[(0)]);
			alg[52] = alg[24]*d[(1)]+alg[32]*(1.0-d[(1)]);
			alg[64] = alg[52];
			alg[68] = alg[48];
			alg[72] = alg[64]*__WSum_3_w[0]+alg[68]*__WSum_3_w[1];
			alg[84] = alg[72];
			alg[88] = 1.0/alg[84];
			alg[92] = alg[88];
			alg[96] = alg[48];
			alg[100] = alg[52];
			alg[104] = alg[88];
			alg[108] = alg[92]*alg[96];
			alg[112] = alg[108];
			alg[116] = alg[52];
			alg[120] = alg[112]*alg[116];
			alg[124] = alg[120];
			alg[128] = x[4];
			alg[132] = alg[124]*alg[128];
			alg[136] = alg[100]*alg[104];
			alg[140] = alg[136];
			alg[144] = alg[132];
			alg[148] = alg[140]*__WSum_11_w[0]+alg[144]*__WSum_11_w[1];
			alg[152] = alg[148];
			alg[156] = x[0];
			alg[184] = alg[152]*__WSum_21_w[0]+alg[156]*__WSum_21_w[1];
			alg[188] = alg[184];
			dx[1] = alg[188];
			return;
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	switch(i)
	{
		case 0:
			alg[0] = __Constant_14_k;
			alg[4] = 100000.0;
			alg[16] = x[0];
			alg[20] = x[4];
			alg[24] = alg[4];
			alg[32] = alg[0];
			alg[36] = alg[0];
			alg[44] = alg[4];
			alg[48] = alg[36]*d[(0)]+alg[44]*(1.0-d[(0)]);
			alg[52] = alg[24]*d[(1)]+alg[32]*(1.0-d[(1)]);
			alg[64] = alg[52];
			alg[68] = alg[48];
			alg[72] = alg[64]*__WSum_3_w[0]+alg[68]*__WSum_3_w[1];
			alg[76] = alg[16]*__WSum_20_w[0]+alg[20]*__WSum_20_w[1];
			alg[80] = alg[76];
			alg[84] = alg[72];
			alg[88] = 1.0/alg[84];
			alg[92] = alg[88];
			alg[96] = alg[48];
			alg[100] = alg[52];
			alg[104] = alg[88];
			alg[108] = alg[92]*alg[96];
			alg[112] = alg[108];
			alg[116] = alg[52];
			alg[120] = alg[112]*alg[116];
			alg[124] = alg[120];
			alg[128] = x[4];
			alg[132] = alg[124]*alg[128];
			alg[136] = alg[100]*alg[104];
			alg[140] = alg[136];
			alg[144] = alg[132];
			alg[148] = alg[140]*__WSum_11_w[0]+alg[144]*__WSum_11_w[1];
			alg[152] = alg[148];
			alg[156] = x[0];
			alg[184] = alg[152]*__WSum_21_w[0]+alg[156]*__WSum_21_w[1];
			alg[188] = alg[184];
			der[0 + 1] = alg[80];
			der[4 + 1] = alg[188];
			return;
		case 1:
			alg[0] = __Constant_14_k;
			alg[4] = 100000.0;
			alg[16] = x[0];
			alg[20] = x[4];
			alg[24] = alg[4];
			alg[32] = alg[0];
			alg[36] = alg[0];
			alg[44] = alg[4];
			alg[48] = alg[36]*d[(0)]+alg[44]*(1.0-d[(0)]);
			alg[52] = alg[24]*d[(1)]+alg[32]*(1.0-d[(1)]);
			alg[64] = alg[52];
			alg[68] = alg[48];
			alg[72] = alg[64]*__WSum_3_w[0]+alg[68]*__WSum_3_w[1];
			alg[76] = alg[16]*__WSum_20_w[0]+alg[20]*__WSum_20_w[1];
			alg[80] = alg[76];
			alg[84] = alg[72];
			alg[88] = 1.0/alg[84];
			alg[92] = alg[88];
			alg[96] = alg[48];
			alg[100] = alg[52];
			alg[104] = alg[88];
			alg[108] = alg[92]*alg[96];
			alg[112] = alg[108];
			alg[116] = alg[52];
			alg[120] = alg[112]*alg[116];
			alg[124] = alg[120];
			alg[128] = x[4];
			alg[132] = alg[124]*alg[128];
			alg[136] = alg[100]*alg[104];
			alg[140] = alg[136];
			alg[144] = alg[132];
			alg[148] = alg[140]*__WSum_11_w[0]+alg[144]*__WSum_11_w[1];
			alg[152] = alg[148];
			alg[156] = x[0];
			alg[184] = alg[152]*__WSum_21_w[0]+alg[156]*__WSum_21_w[1];
			alg[188] = alg[184];
			der[0 + 1] = alg[80];
			der[4 + 1] = alg[188];
			return;
	}
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	switch(i)
	{
		case 0:
			alg[8] = d[(2)];
			alg[40] = alg[8];
			zc[0] = alg[40]-(__qss_switch_12_level);
			return;
		case 1:
			alg[12] = d[(4)];
			alg[28] = alg[12];
			zc[0] = alg[28]-(0.0);
			return;
		case 2:
			zc[0] = t-(d[(3)]);
			return;
		case 3:
			alg[0] = __Constant_14_k;
			alg[4] = 100000.0;
			alg[24] = alg[4];
			alg[32] = alg[0];
			alg[36] = alg[0];
			alg[44] = alg[4];
			alg[48] = alg[36]*d[(0)]+alg[44]*(1.0-d[(0)]);
			alg[52] = alg[24]*d[(1)]+alg[32]*(1.0-d[(1)]);
			alg[56] = alg[52];
			alg[60] = 1.0/alg[56];
			alg[64] = alg[52];
			alg[68] = alg[48];
			alg[72] = alg[64]*__WSum_3_w[0]+alg[68]*__WSum_3_w[1];
			alg[84] = alg[72];
			alg[88] = 1.0/alg[84];
			alg[92] = alg[88];
			alg[96] = alg[48];
			alg[100] = alg[52];
			alg[104] = alg[88];
			alg[108] = alg[92]*alg[96];
			alg[112] = alg[108];
			alg[116] = alg[52];
			alg[120] = alg[112]*alg[116];
			alg[124] = alg[120];
			alg[128] = x[4];
			alg[132] = alg[124]*alg[128];
			alg[136] = alg[100]*alg[104];
			alg[140] = alg[136];
			alg[144] = alg[132];
			alg[148] = alg[140]*__WSum_11_w[0]+alg[144]*__WSum_11_w[1];
			alg[160] = alg[60];
			alg[164] = alg[148];
			alg[168] = alg[160]*alg[164];
			alg[172] = alg[148];
			alg[180] = alg[168];
			alg[192] = alg[172]*d[(5)]+alg[180]*(1.0-d[(5)]);
			alg[196] = alg[192];
			zc[0] = alg[196]-(__hysteretic_18_xu);
			return;
		case 4:
			alg[0] = __Constant_14_k;
			alg[4] = 100000.0;
			alg[24] = alg[4];
			alg[32] = alg[0];
			alg[36] = alg[0];
			alg[44] = alg[4];
			alg[48] = alg[36]*d[(0)]+alg[44]*(1.0-d[(0)]);
			alg[52] = alg[24]*d[(1)]+alg[32]*(1.0-d[(1)]);
			alg[56] = alg[52];
			alg[60] = 1.0/alg[56];
			alg[64] = alg[52];
			alg[68] = alg[48];
			alg[72] = alg[64]*__WSum_3_w[0]+alg[68]*__WSum_3_w[1];
			alg[84] = alg[72];
			alg[88] = 1.0/alg[84];
			alg[92] = alg[88];
			alg[96] = alg[48];
			alg[100] = alg[52];
			alg[104] = alg[88];
			alg[108] = alg[92]*alg[96];
			alg[112] = alg[108];
			alg[116] = alg[52];
			alg[120] = alg[112]*alg[116];
			alg[124] = alg[120];
			alg[128] = x[4];
			alg[132] = alg[124]*alg[128];
			alg[136] = alg[100]*alg[104];
			alg[140] = alg[136];
			alg[144] = alg[132];
			alg[148] = alg[140]*__WSum_11_w[0]+alg[144]*__WSum_11_w[1];
			alg[160] = alg[60];
			alg[164] = alg[148];
			alg[168] = alg[160]*alg[164];
			alg[172] = alg[148];
			alg[180] = alg[168];
			alg[192] = alg[172]*d[(5)]+alg[180]*(1.0-d[(5)]);
			alg[196] = alg[192];
			zc[0] = alg[196]-(__hysteretic_18_xl);
			return;
		case 5:
			alg[12] = d[(4)];
			alg[176] = alg[12];
			zc[0] = alg[176]-(0.0);
			return;
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 0:
			d[(0)] = 1.0;
			return;
		case 1:
			d[(1)] = 1.0;
			return;
		case 2:
			d[(2)] = 1.0-d[(2)];
			d[(3)] = t+d[(2)]*__square_sci_16_DC/10000.0+(1.0-d[(2)])*(1.0-__square_sci_16_DC)/10000.0;
			return;
		case 3:
			d[(4)] = 1.0;
			return;
		case 5:
			d[(5)] = 1.0;
			return;
	}
}

void
MOD_handlerNeg(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 0:
			d[(0)] = 0.0;
			return;
		case 1:
			d[(1)] = 0.0;
			return;
		case 4:
			d[(4)] = __hysteretic_18_yl;
			return;
		case 5:
			d[(5)] = 0.0;
			return;
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	switch(i)
	{
		case 0:
			out[0] = x[4];
			return;
		case 1:
			out[0] = x[0];
			return;
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *discretes = (int*)malloc(6*sizeof(int));
	int *events = (int*)malloc(6*sizeof(int));
	int *outputs = (int*)malloc(2*sizeof(int));
	int *states = (int*)malloc(2*sizeof(int));
	int i26;
	int i27;
	int i28;
	int i29;
	int i;
	simulator->data = QSS_Data(2,6,6,0,50,"buck_disc_qss");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__QSSIntegrator_1_x0 = 0.0;
	__QSSIntegrator_2_x0 = 0.0;
	__qss_switch_12_level = 5.000000000000000000000000e-01;
	__qss_switch_13_level = 0.0;
	__Constant_14_k = 1.000000000000000081803054e-05;
	__Constant_15_k = 100000.0;
	__square_sci_16_amplitude = 1.0;
	__square_sci_16_freq = 10000.0;
	__square_sci_16_DC = 50.0/100.0;
	__hysteretic_18_xl = ((((-9.999999999999999547481118e-07))));
	__hysteretic_18_xu = 9.999999999999999547481118e-07;
	__hysteretic_18_yl = ((((-1.0))));
	__hysteretic_18_yu = 1.0;
	__qss_switch_19_level = 0.0;
	modelData->x[4] = 0.0;
	modelData->x[0] = 0.0;
	modelData->d[(2)] = 1.0;
	modelData->d[(3)] = 0.0;
	// Initialize model code.
	if(modelData->alg[40]>__qss_switch_12_level)
	{
		modelData->d[(0)] = 1.0;
	}
	else if(modelData->alg[40]<__qss_switch_12_level)
	{
		modelData->d[(0)] = 0.0;
	}
	if(modelData->alg[28]>0.0)
	{
		modelData->d[(1)] = 1.0;
	}
	else if(modelData->alg[28]<0.0)
	{
		modelData->d[(1)] = 0.0;
	}
		modelData->d[(3)] = __square_sci_16_DC/10000.0;
	if(modelData->alg[196]>__hysteretic_18_xu)
	{
		modelData->d[(4)] = 1.0;
	}
	if(modelData->alg[196]<__hysteretic_18_xl)
	{
		modelData->d[(4)] = __hysteretic_18_yl;
	}
	if(modelData->alg[176]>0.0)
	{
		modelData->d[(5)] = 1.0;
	}
	else if(modelData->alg[176]<0.0)
	{
		modelData->d[(5)] = 0.0;
	}
		__QSSIntegrator_1_p[(0)] = 0.0;
		__QSSIntegrator_1_p[(1)] = 1.000000000000000055511151e-01;
		__QSSIntegrator_1_p[(2)] = 1.000000000000000020816682e-03;
		__QSSIntegrator_1_p[(3)] = 0.0;
		__QSSIntegrator_2_p[(0)] = 0.0;
		__QSSIntegrator_2_p[(1)] = 1.000000000000000020816682e-03;
		__QSSIntegrator_2_p[(2)] = 1.000000000000000020816682e-03;
		__QSSIntegrator_2_p[(3)] = 0.0;
		__WSum_3_p[(0)] = 1.0;
		__WSum_3_p[(1)] = 1.0;
		__WSum_3_p[(2)] = 0.0;
		__WSum_3_p[(3)] = 0.0;
		__WSum_3_p[(4)] = 0.0;
		__WSum_3_p[(5)] = 0.0;
		__WSum_3_p[(6)] = 0.0;
		__WSum_3_p[(7)] = 0.0;
		__WSum_3_p[(8)] = 2.0;
	for(i26 = 0; i26 <= 1; i26++)
	{
		__WSum_3_w[(i26)] = __WSum_3_p[(i26)];
	}
		__WSum_11_p[(0)] = 24.0;
		__WSum_11_p[(1)] = ((((-1.0))));
		__WSum_11_p[(2)] = 0.0;
		__WSum_11_p[(3)] = 0.0;
		__WSum_11_p[(4)] = 0.0;
		__WSum_11_p[(5)] = 0.0;
		__WSum_11_p[(6)] = 0.0;
		__WSum_11_p[(7)] = 0.0;
		__WSum_11_p[(8)] = 2.0;
	for(i27 = 0; i27 <= 1; i27++)
	{
		__WSum_11_w[(i27)] = __WSum_11_p[(i27)];
	}
		__qss_switch_12_p[(0)] = 5.000000000000000000000000e-01;
		__qss_switch_13_p[(0)] = 0.0;
		__Constant_14_p[(0)] = 1.000000000000000081803054e-05;
		__Constant_15_p[(0)] = 100000.0;
		__square_sci_16_p[(0)] = 1.0;
		__square_sci_16_p[(1)] = 10000.0;
		__square_sci_16_p[(2)] = 50.0;
		__hysteretic_18_p[(0)] = ((((-9.999999999999999547481118e-07))));
		__hysteretic_18_p[(1)] = 9.999999999999999547481118e-07;
		__hysteretic_18_p[(2)] = ((((-1.0))));
		__hysteretic_18_p[(3)] = 1.0;
		__qss_switch_19_p[(0)] = 0.0;
		__WSum_20_p[(0)] = ((((-1000.0))));
		__WSum_20_p[(1)] = 10000.0;
		__WSum_20_p[(2)] = 0.0;
		__WSum_20_p[(3)] = 0.0;
		__WSum_20_p[(4)] = 0.0;
		__WSum_20_p[(5)] = 0.0;
		__WSum_20_p[(6)] = 0.0;
		__WSum_20_p[(7)] = 0.0;
		__WSum_20_p[(8)] = 2.0;
	for(i28 = 0; i28 <= 1; i28++)
	{
		__WSum_20_w[(i28)] = __WSum_20_p[(i28)];
	}
		__WSum_21_p[(0)] = 10000.0;
		__WSum_21_p[(1)] = ((((-10000.0))));
		__WSum_21_p[(2)] = 0.0;
		__WSum_21_p[(3)] = 0.0;
		__WSum_21_p[(4)] = 0.0;
		__WSum_21_p[(5)] = 0.0;
		__WSum_21_p[(6)] = 0.0;
		__WSum_21_p[(7)] = 0.0;
		__WSum_21_p[(8)] = 2.0;
	for(i29 = 0; i29 <= 1; i29++)
	{
		__WSum_21_w[(i29)] = __WSum_21_p[(i29)];
	}
		modelData->d[(5)] = 1.0;
		modelData->d[(4)] = 1.0;
		modelData->d[(2)] = 1.0;
		modelData->d[(1)] = 1.0;
		modelData->d[(0)] = 1.0;
	modelData->nDS[0]++;
	modelData->nDS[0]++;
	modelData->nDS[1]++;
	modelData->nDS[1]++;
	modelData->nSD[0]++;
	modelData->nSD[1]++;
	modelData->nSD[0]++;
	modelData->nSD[1]++;
	modelData->nZS[3]++;
	modelData->nZS[4]++;
	modelData->nSZ[1]++;
	modelData->nSZ[1]++;
	modelData->nHZ[0] += 2;
	modelData->nHZ[1] += 2;
	modelData->nHZ[2] = 1;
	modelData->nHZ[2] += 1;
	modelData->nHZ[3] += 2;
	modelData->nHZ[4] += 2;
	modelData->nHZ[5] += 2;
	modelData->nHD[0] = 1;
	modelData->nHD[1] = 1;
	modelData->event[0].nLHSDsc = 1;
	modelData->event[1].nLHSDsc = 1;
	modelData->event[2].nLHSDsc = 2;
	modelData->event[3].nLHSDsc = 1;
	modelData->event[4].nLHSDsc = 1;
	modelData->event[5].nLHSDsc = 1;
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,2);

	modelData->DS[0][states[0]++] = 0;
	modelData->DS[0][states[0]++] = 1;
	modelData->DS[1][states[1]++] = 0;
	modelData->DS[1][states[1]++] = 1;
		cleanVector(states,0,2);

	modelData->SD[0][states[0]++] = 0;
	modelData->SD[1][states[1]++] = 0;
	modelData->SD[0][states[0]++] = 1;
	modelData->SD[1][states[1]++] = 1;
		cleanVector(events,0,6);

	modelData->ZS[3][events[3]++] = 1;
	modelData->ZS[4][events[4]++] = 1;
		cleanVector(states,0,2);

	modelData->SZ[1][states[1]++] = 3;
	modelData->SZ[1][states[1]++] = 4;
		cleanVector(events,0,6);

	modelData->HZ[0][events[0]++] = 3;
	modelData->HZ[0][events[0]++] = 4;
	modelData->HZ[1][events[1]++] = 3;
	modelData->HZ[1][events[1]++] = 4;
	modelData->HZ[2][events[2]++] = 2;
	modelData->HZ[2][events[2]++] = 0;
	modelData->HZ[3][events[3]++] = 1;
	modelData->HZ[3][events[3]++] = 5;
	modelData->HZ[4][events[4]++] = 1;
	modelData->HZ[4][events[4]++] = 5;
	modelData->HZ[5][events[5]++] = 3;
	modelData->HZ[5][events[5]++] = 4;
		cleanVector(events,0,6);

	modelData->HD[0][events[0]++] = 1;
	modelData->HD[1][events[1]++] = 1;
		cleanVector(events,0,6);

	modelData->event[0].LHSDsc[events[0]++] = 0;
	modelData->event[1].LHSDsc[events[1]++] = 1;
	modelData->event[2].LHSDsc[events[2]++] = 2;
	modelData->event[2].LHSDsc[events[2]++] = 3;
	modelData->event[3].LHSDsc[events[3]++] = 4;
	modelData->event[4].LHSDsc[events[4]++] = 4;
	modelData->event[5].LHSDsc[events[5]++] = 5;
		cleanVector(events,0,6);

	modelData->event[0].direction = 0;
	modelData->event[0].relation = 2;
	modelData->event[1].direction = 0;
	modelData->event[1].relation = 2;
	modelData->event[2].direction = 1;
	modelData->event[2].relation = 2;
	modelData->event[3].direction = 1;
	modelData->event[3].relation = 2;
	modelData->event[4].direction = -1;
	modelData->event[4].relation = 0;
	modelData->event[5].direction = 0;
	modelData->event[5].relation = 2;
	simulator->time = QSS_Time(2,6,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("buck_disc_qss",2,6,2,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

		modelOutput->nOS[0] = 1;
		modelOutput->nSO[1]++;
		modelOutput->nOS[1] = 1;
		modelOutput->nSO[0]++;
	SD_allocOutputMatrix(modelOutput,2,6);
		cleanVector(states,0,2);

		cleanVector(outputs,0,2);

		sprintf(modelOutput->variable[0].name,"QSSIntegrator_1_y[1]");
		sprintf(modelOutput->variable[1].name,"QSSIntegrator_2_y[1]");
		cleanVector(outputs,0,2);

		modelOutput->SO[1][states[1]++] = 0;
		modelOutput->OS[0][outputs[0]++] = 1;
		modelOutput->SO[0][states[0]++] = 1;
		modelOutput->OS[1][outputs[1]++] = 0;
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,MOD_zeroCrossing,MOD_handlerPos,MOD_handlerNeg);
	free(discretes);
	free(events);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
