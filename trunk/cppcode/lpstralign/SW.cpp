#include "SW.h"
#include "math.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

double getSubstituteCost(DATATYPE* a, DATATYPE* b, CModel* model)
{
	double cost = 0;
	for (int i = 0; i < model->m_iMatchCount; i++)
	{
		int t = model->m_vMatch[i];
		int idx = model->m_vFeatureIndex[t];
		int wid = model->m_vWeightIndex[t];
		if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
		{
			cost += model->m_vWeight[wid];
		}
		else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
		{
			DATATYPE fa = a[idx];
			DATATYPE fb = b[idx];
			if (fabs(fa) + fabs(fb) > 0)
			{
				DATATYPE d = fa - fb;
				cost += model->m_vWeight[wid] * d * d / (fabs(fa) + fabs(fb));
				//fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
			}
		}
		else if (model->m_vMapType[t] == MAPTYPE_ABS)
		{
			DATATYPE fa = fabs(a[idx]);
			DATATYPE fb = fabs(b[idx]);
			if (fa + fb > 0)
			{
				DATATYPE d = fa - fb;
				cost += model->m_vWeight[wid] * d * d / (fa + fb);
				//fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
			}
		}
		else if (model->m_vMapType[t] == MAPTYPE_EU)
		{
			DATATYPE fa = a[idx];
			DATATYPE fb = b[idx];
			cost += model->m_vWeight[wid] * fabs(fa - fb);
		}
		else if (model->m_vMapType[t] == MAPTYPE_ABSEU)
		{
			DATATYPE fa = fabs(a[idx]);
			DATATYPE fb = fabs(b[idx]);
			cost += model->m_vWeight[wid] * fabs(fa - fb);
		}
	}
	//double cost = model->m_vWeight[model->m_vMatchMap[model->m_iFeatureDim]] + f;
	//if (cost > 0)
	//fprintf(stderr, "(%f, %f, %d)", cost, model->m_vWeight[model->m_vMatchMap[model->m_iFeatureDim]], model->m_vMatchMap[model->m_iFeatureDim]);
	//fprintf(stderr, "(%f)", cost);
	return cost;
}

double getGapCost(DATATYPE* a, CModel* model)
{
	//fprintf(stderr, "cal gap cost\n");
	double cost = 0;
	for (int i = 0; i < model->m_iGapCount; i++)
	{

		int t = model->m_vGap[i];
		//fprintf(stderr, "t %d\n", t);
		int idx = model->m_vFeatureIndex[t];
		int wid = model->m_vWeightIndex[t];

		if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
		{
			cost += model->m_vWeight[wid];
		}
		else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
		{
			//fprintf(stderr, "gap cost: %d %d %d\n", t, idx, wid);
			DATATYPE fa = a[idx];
			cost += model->m_vWeight[wid] * fabs(fa);
		}
	}
	//double cost = model->m_vWeight[model->m_vGapMap[model->m_iFeatureDim]] + f;
	//fprintf(stderr, "(%f, %f)", cost, model->m_vWeight[model->m_vGapMap[model->m_iFeatureDim]]);
	return cost;
	// return f; //debug
}

double SmithWaterman(CSetOfSeq* pSS1, CSetOfSeq* pSS2, int iSSIndex1,
		int iSSIndex2, CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign,
		CModel* model)
{

	//allocate nodes

	CSWNode ** S;
	//fprintf(stderr, "ttok-1\n");
	int iLengthA = pSeqA->m_iPoint;
	int iLengthB = pSeqB->m_iPoint;
	S = new CSWNode*[iLengthA + 1];
	//fprintf(stderr, "ttok0 length %d %d\n", iLengthA, iLengthB);
	for (int i = 0; i < iLengthA + 1; i++)
	{
		S[i] = new CSWNode[iLengthB + 1];
	}
	for (int i = 0; i < iLengthA + 1; i++)
	{
		S[i][0].score = 0;
		S[i][0].prev = NULL;
		S[i][0].row_index = i;
		S[i][0].column_index = 0;
	}

	for (int i = 0; i < iLengthB + 1; i++)
	{
		S[0][i].score = 0;
		S[0][i].prev = NULL;
		S[0][i].row_index = 0;
		S[0][i].column_index = i;
	}
	//fprintf(stderr, "ttok1\n");
	int maxA = 0, maxB = 0;
	CSWNode* endSequence = &S[0][0];
	DATATYPE** A = pSeqA->m_vFeature;
	DATATYPE** B = pSeqB->m_vFeature;
	//fprintf(stderr, "ttok2\n");
	for (int i = 1; i < iLengthA + 1; i++)
	{
		for (int j = 1; j < iLengthB + 1; j++)
		{
			S[i][j].row_index = i;
			S[i][j].column_index = j;
			S[i][j].score = 0;
			S[i][j].prev = NULL;
			double subst_ij1 = getSubstituteCost(A[i - 1], B[j - 1], model) + S[i - 1][j - 1].score;
			//if (subst_ij1 > 0) fprintf(stderr, "sub cost %f\n", subst_ij1);
			if (subst_ij1 > S[i][j].score)
			{
				S[i][j].score = subst_ij1;
				S[i][j].prev = &S[i - 1][j - 1];
				S[i][j].operation = SUBST;
			}
			double subst_ij2 = getGapCost(A[i - 1], model) + S[i - 1][j].score;
			if (subst_ij2 > S[i][j].score)
			{
				S[i][j].score = subst_ij2;
				S[i][j].prev = &S[i - 1][j];
				S[i][j].operation = DELET;
			}
			double subst_ij3 = getGapCost(B[j - 1], model) + S[i][j - 1].score;
			//fprintf(stderr, "gap cost %f, %f\n", subst_ij2, subst_ij3);
			if (subst_ij3 > S[i][j].score)
			{
				S[i][j].score = subst_ij3;
				S[i][j].prev = &S[i][j - 1];
				S[i][j].operation = INSRT;
			}

			//record information if new maximum
			if (S[i][j].score > endSequence->score)
			{
				endSequence = &S[i][j];
				maxA = i;
				maxB = j;
			}
			//fprintf(stderr, "(%d, %d): %4.2f ", i, j, S[i][j].score);
		}
	}
	//fprintf(stderr, "ttok3\n");
	///////////////////////////////////////////
	// Retrieve alignment
	///////////////////////////////////////////
	/*get length of alignment*/
	CSWNode* temp_node = endSequence;
	pAlign->m_fScore = temp_node->score;
	int oper_count = 0;
	while (temp_node->prev)
	{
		//fprintf(stderr, "oper %d %d\n", oper_count, temp_node->operation);
		pAlign->m_SeqIndex1.push_back(iSSIndex1);
		pAlign->m_SeqIndex2.push_back(iSSIndex2);
		pAlign->m_PointIndex1.push_back(temp_node->prev->row_index);
		pAlign->m_PointIndex2.push_back(temp_node->prev->column_index);
		pAlign->m_operation.push_back(temp_node->operation);
		oper_count++;
		temp_node = temp_node->prev;
	}
	//fprintf(stderr, "ttok4\n");

	/*free DP matrices */
	for (int i = 0; i < iLengthA + 1; i++)
	{
		free(S[i]);
	}
	free(S);
	return pAlign->m_fScore;

}

double SmithWaterman_SetOfSeq(CSetOfSeq* pSS1, CSetOfSeq* pSS2,
		CAlignment* pASet, CModel* model)
{
	CAlignment best;
	//fprintf(stderr, "%d %d\n", pSS1->m_iSeqNum, pSS2->m_iSeqNum);
	for (int i = 0; i < pSS1->m_iSeqNum; i++)
	{
		double maxf = -1;
		int maxind = -1;
		CSequence* pS1 = pSS1->m_vSeqs[i];
		for (int j = 0; j < pSS2->m_iSeqNum; j++)
		{
			CSequence* pS2 = pSS2->m_vSeqs[j];
			CAlignment align;
			align.m_pSS1 = pSS1;
			align.m_pSS2 = pSS2;
			//fprintf(stderr, "ssok 1 %d %d\n", i, j);
			SmithWaterman(pSS1, pSS2, i, j, pS1, pS2, &align, model);
			//fprintf(stderr, "ssok 2\n");
			//fprintf(stderr, "align (%d) to seq (%d), score %f, length %d\n", i, j, align.m_fScore, (int)align.m_operation.size());
			//fprintf(stderr, "\rscore %f, length %d", align.m_fScore, (int)align.m_operation.size());
			if (align.m_fScore > maxf)
			{
				maxf = align.m_fScore;
				maxind = j;
				best = align;
				//fprintf(stderr, "ok j%d\n", j);
			}
		}
		pASet->AddAlignment(best);
		//fprintf(stderr, "ok i%d\n", i);

	}
	//fprintf(stderr, "ok 3\n");

	return pASet->m_fScore;
}

double GetPhi(CAlignment* pAlign, double* phi, int iParamDim, CModel* model)
{
	memset(phi, 0, sizeof(double) * iParamDim);

	int icount = 0;
	for (int i = 0; i < (int) pAlign->m_operation.size(); i++)
	{

		DATATYPE* a = pAlign->m_pSS1->m_vSeqs[pAlign->m_SeqIndex1[i]]->GetPoint(
						pAlign->m_PointIndex1[i]);
		DATATYPE* b = pAlign->m_pSS2->m_vSeqs[pAlign->m_SeqIndex2[i]]->GetPoint(
						pAlign->m_PointIndex2[i]);
		if (pAlign->m_operation[i] == SUBST)
		{
			icount++;
			for (int j = 0; j < model->m_iMatchCount; j++)
			{

				int t = model->m_vMatch[j];
				int idx = model->m_vFeatureIndex[t];
				int wid = model->m_vWeightIndex[t];
				//fprintf(stderr, "match, %d, %d %d \n", t, idx, wid);
				if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
				{
					phi[wid] += 1;
				}
				else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
				{
					DATATYPE fa = a[idx];
					DATATYPE fb = b[idx];
					if (fabs(fa) + fabs(fb) > 0)
					{
						DATATYPE d = fa - fb;
						DATATYPE f = d * d / (fabs(fa) + fabs(fb));
						phi[wid] += f;
						//fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
					}
				}
				else if (model->m_vMapType[t] == MAPTYPE_ABS)
				{
					DATATYPE fa = fabs(a[idx]);
					DATATYPE fb = fabs(b[idx]);
					if (fa + fb > 0)
					{
						DATATYPE d = fa - fb;
						DATATYPE f = d * d / (fa + fb);
						phi[wid] += f;
						//fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
					}
				}
				else if (model->m_vMapType[t] == MAPTYPE_EU)
				{
					DATATYPE fa = a[idx];
					DATATYPE fb = b[idx];
					phi[wid] += fabs(fa - fb);
				}
				else if (model->m_vMapType[t] == MAPTYPE_ABSEU)
				{
					DATATYPE fa = fabs(a[idx]);
					DATATYPE fb = fabs(b[idx]);
					phi[wid] += fabs(fa - fb);
				}

			}
		}
		else if (pAlign->m_operation[i] == DELET || pAlign->m_operation[i] == INSRT)
		{
			icount++;
			for (int j = 0; j < model->m_iGapCount; j++)
			{
				int t = model->m_vGap[j];
				int idx = model->m_vFeatureIndex[t];
				int wid = model->m_vWeightIndex[t];
				//fprintf(stderr, "gap, %d, %d %d \n", t, idx, wid);
				if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
				{
					phi[wid] += 1;
				}
				else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
				{
					DATATYPE fa = a[idx];
					DATATYPE fb = b[idx];
					if (pAlign->m_operation[i] == DELET)
						phi[wid] += fabs(fa);
					else
						phi[wid] += fabs(fb);
				}
			}
		}
		else
		{
			fprintf(stderr, "error operation code %d at %d, %d %d %d %d\n",
					pAlign->m_operation[i], icount, pAlign->m_SeqIndex1[i],
					pAlign->m_SeqIndex2[i], pAlign->m_PointIndex1[i],
					pAlign->m_PointIndex2[i]);
		}
	}
	double fsum = 0;
	for (int i = 0; i < iParamDim; i++)
	{
		fsum += model->m_vWeight[i] * phi[i];
	}
	//fprintf(stderr, "%d %d\n", icount, (int)pAlign->m_operation.size());
	return fsum;
}
