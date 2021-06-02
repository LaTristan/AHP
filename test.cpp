#include<iostream>
#include<cstring>
#include<cstdio>
using namespace std;
struct node{
    int child[256]; //每个节点的儿子数字，child[0]为儿子个数 
    char name[256]; //节点的名称 
    double** comp;  //该节点概率下的下一层成对比较矩阵 
    double CI;      //CI值 
    double* w;      //最大特征值对应的特征向量 
    double p;       //到达该节点的概率 
    double value;   //只对方案层有效，选取该方案的权重 
}* Node;
struct level{
    int num, start_id;//每一层的点数，和起始点的编号 
}* L;
int n_level, n_node; //总层数，总点数 
/*
 * 初始化经验RI数组 
 */
double RI[12]={0,0,0.58,0.9,1.12,1.24,1.32,1.41,1.45,1.49,1.51};
/*
 * 寻找字符c在字符串s中的位置，不存在返回-1
 */ 
int pos(char c, char* s){
    int len = strlen(s);
    for (int i = 0; i < len; i++){
        if (s[i] == c) return i;
    }  
    return -1;
}
/*
 * 读入一个double或者分数形式的数字 
 */
void myScanf(double* p){
    char value[256];
    scanf("%s",value);
    int place;
    if ((place = pos('/', value))>0){
       int f = 0, b = 0;
       for (int i = 0; i < place; i++){
           f = f * 10 + value[i] - '0';
       }
       for (int i = place + 1; i < strlen(value); i++){
           b = b * 10 + value[i] - '0';
       }
       *p = 1.0 * f / b;
    } else{
       *p = 0;
       for (int i = 0; i < strlen(value); i++){
           *p = *p * 10 + value[i] - '0';
       }
    }
}
/*
 * 数据输入过程 
 */
void input()
{
    printf("输入总层数，包括目标层和方案层:\n");
    scanf("%d", &n_level);
    printf("已输入\n");
     
    L = new level[n_level];
    n_node = 0;
    printf("从上往下，输入每层个数:\n");
    for (int i = 0; i < n_level; i++){
        scanf("%d",&L[i].num);
        L[i].start_id = n_node;
        n_node += L[i].num;       
    } 
    printf("已输入\n");
 
    Node = new node[n_node];
    printf("从上往下，输入每层属性名:\n");
    for (int i = 0; i < n_node; i++){
        Node[i].child[0] = 0;
        Node[i].p = 0;
        Node[i].value = 0;
        scanf("%s",Node[i].name);
    }
    printf("已输入\n");
    printf("从上往下，输入层与层之间的关系矩阵:\n");
    for (int i = 0; i < n_level - 1; i++){
        for (int j = 0; j < L[i].num; j++)
        for (int k = 0; k < L[i+1].num; k++){
            int ifConnect;
            scanf("%d",&ifConnect);
            if (ifConnect){
               int cur = L[i].start_id + j;
               Node[cur].child[0]++;
               Node[cur].child[Node[cur].child[0]] = L[i+1].start_id + k;
            }
        }
    }
    printf("已输入\n");
     
    printf("从上往下，输入成对比较矩阵:\n");
    for (int i = 0; i < n_level - 1; i++){
        for (int j = 0; j < L[i].num; j++){
            int cur = L[i].start_id + j;
            int size = Node[cur].child[0];
            Node[cur].comp = new double*[size];
            for (int k = 0; k < size; k++) Node[cur].comp[k] = new double[size];
            for (int li = 0; li < size; li++)
                for (int lj = 0; lj < size; lj++)
                    myScanf(&Node[cur].comp[li][lj]);
        }
    }
    printf("已输入\n");
}
/*
 * 程序中止，并输出原因 
 */
void alert(char* p){
    printf("%s回车退出\n", p);
    freopen("CON", "r", stdin);
    system("pause");
    exit(0);
}
/*
 * n阶矩阵matrix和w相乘，并归一化w
 * retDiv为归一化的倍数，返回归一化后的矩阵 
 */
double* normalize(double** matrix, double* w, int n, double* retDiv){
    double* ret = new double[n];
    for (int i = 0; i < n; i++) ret[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            ret[i] += matrix[i][j] * w[j];
        }
    double sum = 0;
    for (int i = 0; i < n; i++) sum += ret[i];
    for (int i = 0; i < n; i++) ret[i]/=sum;
    *retDiv = sum;
    return ret;
}
/*
 * 检查节点编号为id的节点的成对比较矩阵是否通过一致性检验 
 */
bool checkCR(int id){
    double** matrix = Node[id].comp;
    int n = Node[id].child[0];
    double *w = new double[n], *w_new = new double[n];
    for (int i = 0; i < n; i++) w_new[i] = 1;
    double eps = 1e-10;
    double* retDiv = new double;
    w_new = normalize(matrix, w_new, n, retDiv);
    while (true){
          w = w_new;
          w_new = normalize(matrix, w, n, retDiv);
          bool flag = true;
          for (int i = 0; i < n; i++) 
              if (w[i] - w_new[i] > eps || w[i] - w_new[i] < -eps){
                       flag=false;
                       break;
              }
          if (flag) break;
    }
    Node[id].w = w_new;
    double ans= 0;
    for (int i = 0; i < n; i++)
        ans += w_new[i]/w[i];
    ans = ans * (*retDiv) / n; 
    double CI = (ans - n)/(n-1);
    Node[id].CI = CI;
    if (CI < 0.1 * RI[n]) return true;
    return false;
}
/*
 * 所有单排序权向量的一致性检验 
 */
void checkSingleMatrix(){
    for (int i = 0; i < n_level - 1; i++){
        for (int j = 0; j < L[i].num; j++){
            int cur = L[i].start_id + j;
            if (!checkCR(cur)) {
               printf("第%d层%d号元素<name: %s>的成对比较矩阵无法通过校验,请重新设计",
                   i, j, Node[cur].name);
               alert("");
            }
             
        }
    }
}
/*
 * 总排序权向量的一致性检验 
 */
void checkTotal(){
    double CR = 0;
    Node[0].p = 1;
    /*
     * 计算最后一层中间层的到达概率，即“重要性”
    */
    for (int i = 0; i < n_level - 2; i++){
        for (int j = 0; j < L[i].num; j++){
            int cur = L[i].start_id + j;
            for (int k = 1; k <= Node[cur].child[0]; k++){
                int child = Node[cur].child[k];
                Node[child].p += Node[cur].p * Node[cur].w[k-1];
            }
        }
    }
    int cL = n_level - 2;
    double sumaCI = 0, sumaRI = 0;
    for (int j = 0; j < L[cL].num; j++){
        int cur = L[cL].start_id + j;
        sumaCI += Node[cur].p * Node[cur].CI;
        sumaRI += Node[cur].p * RI[Node[cur].child[0]];
    }
    CR = sumaCI / sumaRI;
    printf("CR值: %.5lf\n",CR);
    if (CR >= 0.1){
        alert("无法通过总排序权向量的一致性检验，请重新设计\n"); 
    }
    printf("检验通过\n"); 
};
/*
 * 根据经检验后的权重，设计最终方案，比较并选择最优解。 
 */
void design(){
    int cL = n_level - 2;
    for (int j = 0; j < L[cL].num; j++){
        int cur = L[cL].start_id + j;
        for (int k = 1; k<= Node[cur].child[0]; k++){
            int child = Node[cur].child[k];
            Node[child].value += Node[cur].p * Node[cur].w[k-1];
        }
    }
    cL = n_level - 1;
    double maxValue = -1e10;
    int maxPlace = -1;
    for (int j = 0; j < L[cL].num; j++){
        int cur = L[cL].start_id + j;
        if (Node[cur].value > maxValue){
           maxValue = Node[cur].value;
           maxPlace = cur;
        }
        printf("<%s>: %.5lf  ", Node[cur].name, Node[cur].value);
    }
    printf("\n\n");
    printf("最优方案为<%s>, 权重为<%.5lf>\n", Node[maxPlace].name, Node[maxPlace].value);
};
/*
 * 整个计算过程 
 */
void calc(){
    printf("单排序权向量一致性检验\n"); 
    checkSingleMatrix();
    printf("检验通过\n"); 
    printf("总排序权向量一致性检验\n"); 
    checkTotal();
    printf("计算方案权值\n");
    design();
    freopen("CON", "r", stdin);
    system("pause");
}
int main()
{
    freopen("level_input.txt","r",stdin);
    input();
    calc();
}
