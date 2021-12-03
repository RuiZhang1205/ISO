#pragma once
#pragma once

#ifndef SHEEP_h
#define SHEEP_h

#include <fstream>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>


#define SHEEP_VERSION "15.1"				//��ǰ�㷨�İ汾��Ϣ������ʱע��ʹ��˫����

using  namespace std;
#define pai acos(-1.0)

#define DIM					30				//����ά��
#define SNUM				40				//��Ⱥ��ģ
#define ITE					300000			//��������
#define N					51				//ÿ��Ĳ�����Ŀ
#define SET					1				//��������
#define DEGREE_HUNTAWAY		0.2				//����Ȯ����̶�
#define Epsilon				1e-20			//GD���
#define step                1E-4			


#define U					0.00000001			//��С��ֵ
#define PRINT_BELLWETHER	1					//�Ƿ����ļ�1�����ÿһ����ͷ��
#define PRINT_HUNTAWAY		1					//�Ƿ��������Ȯ����
#define TIMES_HUNTAWAY		0	 				//����Ȯ���������С���
#define PRINTF_RESULT_FILE	1					//�Ƿ�����ļ�2 �����ÿ�����ս��

/*������Ӧ�Ⱥ���ѡ��Ҫ���ĸ���������Ϊ1,��ֻ����һ��Ϊ1*/
#define TEST_FUN_CHOICE		"shifted_rotated_BentCigar"		//�ı���Ժ�������ͬʱ���˺궨���Ϊ��Ӧ���Ժ�����

#define Best_fitness		100.00						//�ı���Ժ�������ͬʱ���Ĵ˺�����ȫ�����Ž�

#define left_range			-100.00						//ÿһά�����귶Χ(����������������дΪ**.0����������600дΪ600.0)
#define right_range		     100.00						 //ÿһά�����귶Χ(����������������дΪ**.0����������600дΪ600.0)


//CEC2017����
#define shifted_rotated_BentCigar					1	//range=100,bestfitness=100
#define shifted_rotated_sum_of_Different_Power		0	//range=100,bestfitness=200
#define shifted_rotated_Zakharov					0	//range=100,bestfitness=300
#define shifted_rotated_Rosenbrock					0	//range=100,bestfitness=400
#define shifted_rotated_Rastrigin					0	//range=100,bestfitness=500
#define shifted_Rotated_Scaffer7					0	//range=100,bestfitness=600
#define shifted_Rotated_Lunacek_Bi_Rastrigin		0	//range=100,bestfitness=700
#define shifted_Rotated_Non_Continuous_Rastrigin	0	//range=100,bestfitness=800
#define Shifted_Rotated_Levy						0	//range=100,bestfitness=900
#define Shifted_Rotated_Schwefel					0	//range=100,bestfitness=1000

#define hybrid_function_1							0	//range=100,bestfitness=1100
#define hybrid_function_2							0	//range=100,bestfitness=1200
#define hybrid_function_3							0	//range=100,bestfitness=1300
#define hybrid_function_4							0	//range=100,bestfitness=1400
#define hybrid_function_5							0	//range=100,bestfitness=1500
#define hybrid_function_6							0	//range=100,bestfitness=1600
#define hybrid_function_7							0	//range=100,bestfitness=1700
#define hybrid_function_8							0	//range=100,bestfitness=1800
#define hybrid_function_9							0	//range=100,bestfitness=1900
#define hybrid_function_10							0	//range=100,bestfitness=2000

#define Composition_function_1						0	//range=100,Best_fitness=2100
#define Composition_function_2						0	//range=100,Best_fitness=2200
#define Composition_function_3						0	//range=100,Best_fitness=2300
#define Composition_function_4						0	//range=100,Best_fitness=2400
#define Composition_function_5						0	//range=100,Best_fitness=2500
#define Composition_function_6						0	//range=100,Best_fitness=2600
#define Composition_function_7						0	//range=100,Best_fitness=2700
#define Composition_function_8						0	//range=100,Best_fitness=2800
#define Composition_function_9						0	//range=100,Best_fitness=2900
#define Composition_function_10						0	//range=100,Best_fitness=3000


class SHEEP
{
public:
	double	coordinate[DIM];					//�洢ÿֻ�������
	double	fitness;							//�洢�ʶ�ֵ
	int		number;								//�洢��ı��
	int		scatter;							//�ж���ֻ���Ƿ񱻴�ɢ 1����ɢ��0û����ɢ
public:
	friend class GROUPSHEEP;
};
class GROUPSHEEP
{
private:										//������Χ��Ҳ���ǽ�ռ䷶Χ
	double	left;								//�洢ÿά����ķ�Χ
	double	right;
public:
	SHEEP	sheep[SNUM];						//��Ⱥ
	int		bellwethernumber;					//��ͷ����
	double	worstfitness;						//��Ⱥ����ʶ�ֵ
	double	meanfitness;						//������Ⱥ���ʶ�ֵ��ƽ��ֵ
	double	oldbellwetherfitness;				//�洢��һ����ͷ����ʶ�ֵ
	int		generationTimes;					//��ǰ��������
	double  GBest[DIM];							//�������λ��
	double  BeforeGD;							//�ݶ��½�֮ǰ��ֵ
	double  AfterGD;							//�ݶ��½����ֵ
public:
	GROUPSHEEP();
	void	initofgroup();						//��ʼ����Ⱥ
	void	leader();							//��ͷ��׶�
	void	wander();							//��Ⱥ���ν׶�
	void	bellwether();						//����һ����ͷ��,���������Ⱥ��ƽ���ʶ�ֵ������ʶ�ֵ
	int		huntaway();							//����Ȯ�׶�
	void	GradientDescent();					//�ݶ��½�
};

double Computafitness(double a[]);				//���Ժ���		
#endif
#pragma once

