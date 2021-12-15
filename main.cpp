#pragma warning(disable:4996)
#include <iostream>
#include "sheep.h"

FILE *fp;
FILE *fp2;

int main() {
	struct tm *newtime;
	char tmpbuf[128];
	char tmpbuf2[128];
	time_t lt1;
	lt1 = time(NULL);
	newtime = localtime(&lt1);
	srand((unsigned)time(NULL));

	sprintf(tmpbuf, "SHEEPv%s_U%-.2EDIM%d_SNUM%d_ITE%d_TH%d_%d%d_%s_%d%02d%02d_%02d%02d%02d_File1.txt", SHEEP_VERSION, U,
		DIM, SNUM, ITE, TIMES_HUNTAWAY, PRINT_BELLWETHER, PRINT_HUNTAWAY, TEST_FUN_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);
	sprintf(tmpbuf2, "SHEEPv%s_U%-.2EDIM%d_SNUM%d_ITE%d_TH%d_%d%d_%s_%d%02d%02d_%d%02d%02d_File2.txt", SHEEP_VERSION, U,
		DIM, SNUM, ITE, TIMES_HUNTAWAY, PRINT_BELLWETHER, PRINT_HUNTAWAY, TEST_FUN_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);

#if PRINT_BELLWETHER
	fp = fopen(tmpbuf, "wt");
	fprintf(fp, "%s\n\n", tmpbuf);
	/*fprintf(fp, "����ά��:\t%d\n��Ⱥ��ģ:\t%d\n��������:\t%d\n�����귶Χ:\t%f\n�����귶Χ:\t%f\nÿ�������Ŀ:\t%d\n��������:\t%d\n����Ȯ����̶�:\t%.2f\n����Ȯ������С���:\t%d\n����Ȯ���뷧ֵ:\t%lf\n���Ժ�������:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);*/
	fprintf(fp, "����ά��:\t%d\n��Ⱥ��ģ:\t%d\n��������:\t%d\n�����귶Χ:\t%f\n�����귶Χ:\t%f\nÿ�������Ŀ:\t%d\n��������:\t%d\n����Ȯ������С���:\t%d\n����Ȯ���뷧ֵ:\t%lf\n���Ժ�������:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);
#endif
#if PRINTF_RESULT_FILE
	fp2 = fopen(tmpbuf2, "wt");
	fprintf(fp2, "%s\n\n", tmpbuf2);
	/*fprintf(fp2, "����ά��:\t%d\n��Ⱥ��ģ:\t%d\n��������:\t%d\n�����귶Χ:\t%f\n�����귶Χ:\t%f\nÿ�������Ŀ:\t%d\n��������:\t%d\n����Ȯ����̶�:\t%.2f\n����Ȯ������С���:\t%d\n����Ȯ���뷧ֵ:\t%lf\n���Ժ�������:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, DEGREE_HUNTAWAY, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);//���ļ���������β��Եı�ע��Ϣ*/
	fprintf(fp2, "����ά��:\t%d\n��Ⱥ��ģ:\t%d\n��������:\t%d\n�����귶Χ:\t%f\n�����귶Χ:\t%f\nÿ�������Ŀ:\t%d\n��������:\t%d\n����Ȯ������С���:\t%d\n����Ȯ���뷧ֵ:\t%lf\n���Ժ�������:\t%s\n\n\n",
		DIM, SNUM, ITE, left_range, right_range, N, SET, TIMES_HUNTAWAY, U, TEST_FUN_CHOICE);//���ļ���������β��Եı�ע��Ϣ
#endif
	for (int i = 0; i < SET; i++) {
#if PRINT_BELLWETHER
		fprintf(fp, "=======================================================================================================\n");
		fprintf(fp, "��%d��\n", i + 1);
#endif
#if PRINTF_RESULT_FILE
		fprintf(fp2, "======================================================================================================\n");
		fprintf(fp2, "��%d��\n", i + 1);
		fprintf(fp2, "���Դ���\t�ﵽ�������ֵʱ�Ĵ���\t����ֵ\t�������\t���ձ��\t���պ���ֵ\t���պ������\n");
		fprintf(fp2, "ע���ﵽ�������ֵʱ�Ĵ��� ��������Ϊ%5dʱΪδ�������������������ҵ���������ֵ��\n", ITE);
#endif
		long iteSUM = 0, bestITE = LONG_MAX, worstITE = LONG_MIN;
		clock_t start, finish;			//�����������ʱ��
		double sum = 0, best = LONG_MAX, worst = LONG_MIN, duration, satisfactFitness;
		double satisfactSum = 0, satisfactBest = LONG_MAX, satisfactWorst = LONG_MIN;
		start = clock();
		int k = 0;
		int isHuntaway, satisfactGnrtTm, satisfactTms = 0;
		for (k = 0; k < N; k++) {		//��������N �γ��� �õ�N������
			GROUPSHEEP group;
			isHuntaway = 0;
			satisfactGnrtTm = 0;
			group.initofgroup();		//��ʼ������Ⱥ
			group.bellwether();			//�ҳ���ͷ��
#if PRINT_BELLWETHER
			fprintf(fp, "����\t���\t����ֵ\t�������\t�Ƿ�����Ȯ����\t����\t\n");
			fprintf(fp, "%6d\t%3d\t%-E\t%-E\t", group.generationTimes, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
			fprintf(fp, "        \t");
			fprintf(fp, "\n");
#endif
			while (group.generationTimes < ITE) {

				group.leader();					//��ͷ��׶�
				group.wander();					//��Ⱥ����
				group.oldbellwetherfitness = group.sheep[group.bellwethernumber - 1].fitness;
				group.bellwether();				//�ҳ���ͷ��
				isHuntaway = group.huntaway();	//����Ȯ�׶ο�ʼ
				group.bellwether();				//���¸�����ͷ��
												//printf("����Ȯ�׶���ͷ�򣺱��%d    %.6lf\n", group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness);
												//printf("worstfitness:%.6lf\nmeanfitness:%.6lf\n\n\n", group.worstfitness, group.meanfitness);
#if PRINT_BELLWETHER
												//fprintf(fp, "%6d\t%3d\t%-E\t%-E\t", group.generationTimes, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
				fprintf(fp, "%6d\t%3d\t%-E\t%-E\t", group.generationTimes, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#if PRINT_HUNTAWAY
				if (isHuntaway == 1)
				{
					//fprintf(fp, "Huntaway\tbefore GD:%E\t\tafter GD:%E\t", group.BeforeGD, group.AfterGD);
					fprintf(fp, "Huntaway");
				}
#endif
				if (isHuntaway == 1)
					fprintf(fp, "\t");
				else
					fprintf(fp, "        \t");
				//for (int dim = 0;dim < DIM;dim++)
				//fprintf(fp, "%-E\t", group.sheep[group.bellwethernumber - 1].coordinate[dim]);
				fprintf(fp, "\n");
#endif
				if (group.sheep[group.bellwethernumber - 1].fitness - Best_fitness < 1E-06 && satisfactGnrtTm == 0 || group.generationTimes == (ITE) && satisfactGnrtTm == 0) {
					iteSUM = iteSUM + (group.generationTimes);											//��¼���������ĺ�
					if ((group.generationTimes) < bestITE) bestITE = group.generationTimes;				//��¼���ŵ�������
					if ((group.generationTimes) > worstITE) worstITE = group.generationTimes;			//��¼���ŵ�������
					satisfactGnrtTm = group.generationTimes;
					satisfactFitness = group.sheep[group.bellwethernumber - 1].fitness;
					if (satisfactGnrtTm < ITE) satisfactTms++;
				}
			}
			satisfactSum = satisfactFitness;			//����20�����ݵĺ�
			if (satisfactFitness < satisfactBest)			//�ҳ�N�����������õ�ֵ
				satisfactBest = satisfactFitness;
			if (satisfactFitness > satisfactWorst)		//�ҳ�N������������ֵ
				satisfactWorst = satisfactFitness;

			sum = sum + group.sheep[group.bellwethernumber - 1].fitness;	//����20�����ݵĺ�
			if (group.sheep[group.bellwethernumber - 1].fitness < best)		//�ҳ�N�����������õ�ֵ
				best = group.sheep[group.bellwethernumber - 1].fitness;
			if (group.sheep[group.bellwethernumber - 1].fitness > worst)		//�ҳ�N������������ֵ
				worst = group.sheep[group.bellwethernumber - 1].fitness;
#if PRINT_BELLWETHER
			if (satisfactGnrtTm == ITE) {
				fprintf(fp, "��%d�β��ԣ�\n��\t%5d\t��ֹͣδ�����������������ҵ��������ĺ���ֵ��", k + 1, satisfactGnrtTm);
			}
			else
				fprintf(fp, "��%d�β��ԣ�\n��\t%5d\t������ֵ�ﵽ�������", k + 1, satisfactGnrtTm);
			fprintf(fp, "����ֵ��\t%E\t������\t%E\n", satisfactFitness, satisfactFitness - Best_fitness);
			fprintf(fp, "��\t%5d\t���ﵽ��ֹ������", group.generationTimes);
			fprintf(fp, "��ͷ��ı��\t%3d\t����ֵ��\t%E\t������\t%E\n", group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#endif
#if PRINTF_RESULT_FILE
			fprintf(fp2, "%d\t%5d\t%E\t%E\t", k + 1, satisfactGnrtTm, satisfactFitness, satisfactFitness - Best_fitness);
			fprintf(fp2, "%3d\t%E\t%E\n", group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness, group.sheep[group.bellwethernumber - 1].fitness - Best_fitness);
#endif
			//printf("%d ���%d    %.6lf\n", k + 1, group.bellwethernumber, group.sheep[group.bellwethernumber - 1].fitness);
			printf("Set %d/%d %d'/%d' is OK.\n", i + 1, SET, k + 1, N);
		}
#if PRINT_BELLWETHER
		if (satisfactTms == 0) {
			fprintf(fp, "\nδ�������޵����������ҵ��������Χ�ڵ�ֵ��\n");
		}
		else {
			fprintf(fp, "\n��\t%d\t�β����У���\t%d\t�β����ҵ�����������ֵ��������������ۺ�ͳ�����£�\n", N, satisfactTms);
			fprintf(fp, "���ź���ֵ��������\t%E\t%E\n", satisfactBest, satisfactBest - Best_fitness);							//N����������������ֵ�����ĵ�
			fprintf(fp, "����ֵ��������\t%E\t%E\n", satisfactWorst, satisfactWorst - Best_fitness);						//N����������������ֵ�����ĵ�
			fprintf(fp, "ƽ������ֵ��������\t%E\t%E\n", satisfactSum / N, satisfactSum / N - Best_fitness);					//N�������������ƽ��ֵ�����ĵ�
			fprintf(fp, "��С����������\t%6d\n", bestITE);				//N���������������õ������������ĵ�
			fprintf(fp, "������������\t%6d\n", worstITE);				//N����������������������������ĵ�
			fprintf(fp, "ƽ������������\t%.1lf\n", (double)iteSUM / k);	//N�������������ƽ���������������ĵ�
		}
#endif
#if PRINTF_RESULT_FILE
		if (satisfactTms == 0) {
			fprintf(fp2, "\nδ�������޵����������ҵ��������Χ�ڵ�ֵ��\n");
		}
		else {
			fprintf(fp2, "\n��\t%d\t�β����У���\t%d\t�β����ҵ�����������ֵ��������������ۺ�ͳ�����£�\n", N, satisfactTms);
			fprintf(fp2, "���ź���ֵ��������\t%E\t%E\n", satisfactBest, satisfactBest - Best_fitness);						//N����������������ֵ�����ĵ�
			fprintf(fp2, "����ֵ��������\t%E\t%E\n", satisfactWorst, satisfactWorst - Best_fitness);						//N����������������ֵ�����ĵ�
			fprintf(fp2, "ƽ������ֵ��������\t%E\t%E\n", satisfactSum / N, satisfactSum / N - Best_fitness);					//N�������������ƽ��ֵ�����ĵ�
			fprintf(fp2, "��С����������\t%6d\n", bestITE);				//N���������������õ������������ĵ�
			fprintf(fp2, "������������\t%6d\n", worstITE);				//N����������������������������ĵ�
			fprintf(fp2, "ƽ������������\t%.1lf\n", (double)iteSUM / k);	//N�������������ƽ���������������ĵ�
		}
#endif
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
#if PRINT_BELLWETHER
		fprintf(fp, "\nÿ�β��Զ�������������������ͳ���������£�\n");
		fprintf(fp, "���ź���ֵ��������\t%E\t%E\n", best, best - Best_fitness);						//N����������������ֵ�����ĵ�
		fprintf(fp, "����ֵ��������\t%E\t%E\n", worst, worst - Best_fitness);						//N����������������ֵ�����ĵ�
		fprintf(fp, "ƽ������ֵ��������\t%E\t%E\n", sum / N, sum / N - Best_fitness);					//N�������������ƽ��ֵ�����ĵ�
		fprintf(fp, "������ʱ�䣺\t%lf��\n", duration);
		fprintf(fp, "======================================================================================================\n\n\n\n");
#endif
#if PRINTF_RESULT_FILE
		fprintf(fp2, "\nÿ�β��Զ�������������������ͳ���������£�\n");
		fprintf(fp2, "���ź���ֵ��������\t%E\t%E\n", best, best - Best_fitness);						//N����������������ֵ�����ĵ�
		fprintf(fp2, "����ֵ��������\t%E\t%E\n", worst, worst - Best_fitness);					//N����������������ֵ�����ĵ�
		fprintf(fp2, "ƽ������ֵ��������\t%E\t%E\n", sum / N, sum / N - Best_fitness);				//N�������������ƽ��ֵ�����ĵ�
		fprintf(fp2, "������ʱ�䣺\t%lf��\n", duration);
		fprintf(fp2, "======================================================================================================\n\n\n\n");
#endif
		//printf("best :%.6lf\n", best);							//N����������������ֵ
		//printf("worst :%.6lf\n", worst);							//N����������������ֵ
		printf("ƽ����Ӧ��: %E\n", sum / N - Best_fitness);			//N�������������ƽ��ֵ
		printf("Set %d is OK.\n", i + 1);
	}
	fclose(fp);
	fclose(fp2);
	cout << "Finish" << endl;
	return 0;
}
