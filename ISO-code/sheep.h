#pragma once
#pragma once

#ifndef SHEEP_h
#define SHEEP_h

#include <fstream>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <cmath>


#define SHEEP_VERSION "15.1"				//当前算法的版本信息，更改时注意使用双引号

using  namespace std;
#define pai acos(-1.0)

#define DIM					30				//粒子维度
#define SNUM				40				//种群规模
#define ITE					300000			//迭代次数
#define N					51				//每组的测试数目
#define SET					1				//测试组数
#define DEGREE_HUNTAWAY		0.2				//牧羊犬介入程度
#define Epsilon				1e-20			//GD误差
#define step                1E-4			


#define U					0.00000001			//最小阀值
#define PRINT_BELLWETHER	1					//是否在文件1中输出每一代领头羊
#define PRINT_HUNTAWAY		1					//是否输出牧羊犬介入
#define TIMES_HUNTAWAY		0	 				//牧羊犬介入代数最小间隔
#define PRINTF_RESULT_FILE	1					//是否输出文件2 仅输出每次最终结果

/*各种适应度函数选择，要用哪个，就设置为1,但只能有一个为1*/
#define TEST_FUN_CHOICE		"shifted_rotated_BentCigar"		//改变测试函数，请同时将此宏定义改为对应测试函数名

#define Best_fitness		100.00						//改变测试函数，请同时更改此函数的全局最优解

#define left_range			-100.00						//每一维度坐标范围(如果是整数，请务必写为**.0，比如整数600写为600.0)
#define right_range		     100.00						 //每一维度坐标范围(如果是整数，请务必写为**.0，比如整数600写为600.0)


//CEC2017函数
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
	double	coordinate[DIM];					//存储每只羊的坐标
	double	fitness;							//存储适度值
	int		number;								//存储羊的编号
	int		scatter;							//判断这只羊是否被打散 1被打散，0没被打散
public:
	friend class GROUPSHEEP;
};
class GROUPSHEEP
{
private:										//牧场范围，也就是解空间范围
	double	left;								//存储每维坐标的范围
	double	right;
public:
	SHEEP	sheep[SNUM];						//羊群
	int		bellwethernumber;					//领头羊编号
	double	worstfitness;						//种群最差适度值
	double	meanfitness;						//整个羊群的适度值的平均值
	double	oldbellwetherfitness;				//存储上一代领头羊的适度值
	int		generationTimes;					//当前迭代次数
	double  GBest[DIM];							//最优羊的位置
	double  BeforeGD;							//梯度下降之前的值
	double  AfterGD;							//梯度下降后的值
public:
	GROUPSHEEP();
	void	initofgroup();						//初始化种群
	void	leader();							//领头羊阶段
	void	wander();							//羊群漫游阶段
	void	bellwether();						//更新一次领头羊,并且求出羊群的平均适度值和最差适度值
	int		huntaway();							//牧羊犬阶段
	void	GradientDescent();					//梯度下降
};

double Computafitness(double a[]);				//测试函数		
#endif
#pragma once

