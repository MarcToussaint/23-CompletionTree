#include <Search/ComputeNode.h>
#include <Search/AStar.h>

#include <Core/array.h>
#include <Core/util.h>
#include <Core/graph.h>
#include <math.h>

#include <LGP/LGP_Tool.h>

//===========================================================================

struct Problem_Options {
  enum ProblemType { none=0, Random, SingularLine, Bandits } type;
  RAI_PARAM("Tree/", int, maxDepth, 3)
  RAI_PARAM("Tree/", int, hiddenTarget, 0)
  RAI_PARAM("Tree/", int, branching, 3)
  RAI_PARAM("Tree/", double, badC, 10)
  RAI_PARAM("Tree/", int, problem, Random)

  RAI_PARAM("ELS/", double, w0, 1.)
  RAI_PARAM("ELS/", double, wP, 2.)
  RAI_PARAM("ELS/", double, c0, 1.)
  RAI_PARAM("ELS/", double, cP, 2.)

  RAI_PARAM("", rai::String, solver, "ELS")

  std::shared_ptr<ofstream> fil;

  Problem_Options() { type = (ProblemType)problem; }
};

//===========================================================================

namespace rai {
  struct ComputeNode_Backup : ComputeNode {
    double y_now=0.;
    double y_tot=0.;
    double y_num=0.;
    double y_best=-1.;

    std::shared_ptr<Problem_Options> opt;

    ComputeNode_Backup(ComputeNode_Backup* parent) : ComputeNode(parent) {}

    virtual double treePolicyScore(int i){
      double beta=1.;
      double score=0.;
      ComputeNode_Backup* par = dynamic_cast<ComputeNode_Backup*>(parent);
      if(y_num){ //sampled before -> normal UCB score
        if(!isTerminal){
          score = y_tot/y_num + beta * sqrt(2. * ::log(par->y_num)/y_num);
        }else{
          score = -1.; //terminals must only be sampled once
        }
      }else{ //never sampled
        if(isComplete){
          score = 1e6; // available but never sampled -> UCB rule
        }else{ //incomplete
          double eff_this = ::pow(c/opt->c0, opt->cP) + ::pow(double(i)/opt->w0, opt->wP);
          double eff_parent = ::pow((par->c_tot+par->y_num)/opt->c0, opt->cP) + 0.; // + ::pow(double(i)/opt->w0, opt->wP);
          if(eff_this < eff_parent) score = 1e6 - eff_this;
          else score = 0.;
        }
      }
      return score;
    }


    void backup_y(double y){
      ComputeNode_Backup *n = this;
      while(n){
        n->y_now = y;
        n->y_tot += y;
        n->y_num++;
        if(n->y_best<0. || y>n->y_best) n->y_best=y;
        n = dynamic_cast<ComputeNode_Backup*>(n->parent);
      }
    }

    virtual void data(Graph& g) const {
      ComputeNode::data(g);
      g.add("C", c_tot);
      g.add("y_n", y_num);
      if(y_num) g.add("<y>",y_tot/y_num);
    }
  };
}

//===========================================================================

//-- a random tree problem
struct HardBandit : rai::ComputeNode_Backup {
  double hat_c=-1.;
  double hat_y=-1.;
  int depth=0;
  bool goodLine=true;


  HardBandit(HardBandit *parent, uint decision) : ComputeNode_Backup(parent) {
    if(!parent){
      opt=make_shared<Problem_Options>();
      if(opt->solver=="RoundRobin"){ opt->c0=opt->cP=opt->w0=opt->wP=1.; }
      depth = 0;
    }else{
      opt = parent->opt;
      depth = parent->depth+1;
    }

    isTerminal = (depth >= opt->maxDepth);

    if(parent){
      goodLine = ((int)decision==opt->hiddenTarget) && dynamic_cast<HardBandit*>(parent)->goodLine;
    }else{ //root
      goodLine = true;
      hat_c=hat_y=l=0.;
      isComplete=true;
    }

    if(parent){
      if(opt->type==opt->Random || opt->type==opt->Bandits){
        hat_c = rnd.uni(1., 10.);
        hat_y = rnd.uni(0., 1.);
        if(parent) hat_y += dynamic_cast<HardBandit*>(parent)->hat_y; //SUM of all hat_y along path
        if(opt->type==opt->Bandits && isTerminal) hat_c=.1; //Bandit samples have fixed compute effort
      }else if(opt->type==opt->SingularLine){
        hat_c = (goodLine? opt->badC : 1.);
        hat_y = (goodLine?1.:0.);
      }else NIY;
    }

    name <<(goodLine?'+':'-') <<'d' <<depth <<'i' <<decision;
  }

  void compute(){
    CHECK(!isComplete, "is already complete!");
    if(c+1.>=hat_c){
      c_now = hat_c - c;
      c = hat_c;
      isComplete=true;
      l=0.;
      if(isTerminal){
        double y = hat_y;
        if(opt->type==opt->Random || opt->type==opt->Bandits) y /= double(depth);
        backup_y(y);
      }
    }else{
      c_now = 1.;
      c += c_now;
    }
    backup_c(c_now);
    f_prio = baseLevel + computePenalty();
  }

  virtual double computePenalty(){
    return ::pow(c/opt->c0, opt->cP);
  }

  double branchingPenalty_child(int i){
    return ::pow(double(i)/opt->w0, opt->wP);
  }

  virtual int getNumDecisions(){
    if(opt->type==opt->Bandits){
      if(depth==opt->maxDepth-1) return -1;
      else return opt->branching;
    }else{
      return opt->branching;
    }
  }

  virtual shared_ptr<ComputeNode> createNewChild(int i){
    CHECK(isComplete, "can't create child of non-complete");
    CHECK(!isTerminal, "can't create child of terminal");
    return make_shared<HardBandit>(this, i);
  }

};

//===========================================================================

void procEval(const char* prefix, uint K){
  arrA X(K);
  for(uint k=0;k<K;k++){
    FILE(STRING(prefix <<k <<".dat")) >>X(k);
  }

  rai::String name=prefix; name <<"VAR.dat";
  ofstream out(name);

  for(uint t=0;t<X(0).d0;t++){
    for(uint i=0;i<X(0).d1;i++){
      double m=0., v=0.;
      for(uint k=0;k<K;k++){
        double x = X(k)(t,i);
        m += x;
        v += x*x;
      }
      m /= double(K);
      v = ::sqrt(v/double(K) - m*m + 1e-10);
      v /= ::sqrt(double(K));
      out <<m <<' ' <<v <<' ';
    }
    out <<endl;
  }
}

//===========================================================================

void ex3_bandit(){
  Problem_Options opt;
  CHECK_EQ(opt.problem, 3, "");

  uint evalLimit = rai::getParameter<double>("evalLimit");

  uint K=10;
  for(uint k=0;k<K;k++){
    std::shared_ptr<HardBandit> root = std::make_shared<HardBandit>(nullptr,0);
    rai::AStar astar(root);
    if(root->opt->solver=="TreePolicy") astar.mode=astar.treePolicy;
    else if(root->opt->solver=="RoundRobin") astar.mode=astar.FIFO;

    ofstream fil(STRING("ex3_data/" <<k <<".dat"));
    uint c_tot=0;
    while(root->c_tot<evalLimit && astar.queue.N){
      astar.step();
      if(floor(root->c_tot)>c_tot){ //write into the file only when a new c_tot integer is reached -- for comparability
        c_tot = floor(root->c_tot);
        // step  c_total y_now y_avg y_best n_nodes n_sol n_level1
        fil <<astar.steps <<' ' <<root->c_tot <<' ' <<root->y_now <<' ' <<(root->y_num>0.?root->y_tot/root->y_num:0.) <<' ' <<root->y_best <<' ' <<astar.mem.N <<' ' <<astar.solutions.N <<' ' <<root->children.N <<endl; //<<' ' <<groundTruthMean <<' ' <<groundTruthBest <<' ' <<terminals.N <<' ' <<nonTerminals.N <<endl;
      }
      //    printTree(astar.mem); rai::wait();
    }
  }
  procEval("ex3_data/", K);
  if(opt.solver=="ELS") opt.solver <<opt.wP <<opt.cP;
  if(opt.solver=="TreePolicy") opt.solver <<opt.wP <<opt.cP;
  rai::system(STRING("cp ex3_data/VAR.dat ex3_data/Bandit_" <<opt.solver <<".dat"));
}

//===========================================================================

void ex2_Random(){
  Problem_Options opt;
  CHECK_EQ(opt.problem, 1, "");
  uint K=10;
  uint evalLimit = rai::getParameter<double>("evalLimit");
  for(uint k=0;k<K;k++){
    auto root = std::make_shared<HardBandit>(nullptr,0);
    rai::AStar astar(root);
    if(root->opt->solver=="TreePolicy") astar.mode=astar.treePolicy;
    else if(root->opt->solver=="RoundRobin") astar.mode=astar.FIFO;

    //ComputeTree_Solver S(root);
    astar.verbose = 0;
    ofstream fil(STRING("ex2_data/" <<k <<".dat"));
    uint c_tot=0;
    while(root->c_tot<evalLimit && astar.queue.N){
      astar.step();
      if(floor(root->c_tot)>c_tot){ //write into the file only when a new c_tot integer is reached -- for comparability
        c_tot = floor(root->c_tot);
        // step  c_total y_now y_avg y_best n_nodes
        fil <<astar.steps <<' ' <<root->c_tot <<' ' <<root->y_now <<' ' <<(root->y_num>0.?root->y_tot/root->y_num:0.) <<' ' <<root->y_best <<' ' <<astar.mem.N <<endl; //<<' ' <<groundTruthMean <<' ' <<groundTruthBest <<' ' <<terminals.N <<' ' <<nonTerminals.N <<endl;
      }
//      printTree(astar.mem); rai::wait();
    }
  }
  procEval("ex2_data/", K);
  gnuplot("load 'plt'");
  if(opt.solver=="ELS") opt.solver <<opt.wP <<opt.cP;
  rai::system(STRING("cp ex2_data/VAR.dat ex2_data/Random_" <<opt.solver <<".dat"));
  rai::wait();
}

//===========================================================================

void evalTimeToTarget(uint K, double& m, double& v){
  uint evalLimit = rai::getParameter<double>("evalLimit");
  m=0., v=0.;
  for(uint k=0;k<K;k++){
    auto root = make_shared<HardBandit>(nullptr,0);
    rai::AStar astar(root);
    if(root->opt->solver=="TreePolicy") astar.mode=astar.treePolicy;
    else if(root->opt->solver=="RoundRobin") astar.mode=astar.FIFO;

    astar.verbose = 0;
    root->opt->fil = make_shared<ofstream>("VAR.dat");
    while(root->y_now<1. && root->c_tot<evalLimit){
      astar.step();
      //printTree(astar.mem);      rai::wait();
    }
    m += root->c_tot;
    v += root->c_tot*root->c_tot;
    //LOG(0) <<"Total cost to target: " <<root->c_tot;
  }

//  printTree(cout, S.root);
//  gnuplot("plot 'VAR.dat' us 2:3 t 'y' w p, '' us 2:4 t 'y_{mean}' lw 2, '' us 2:7 t 'y_{max}' lw 2, '' us 2:(log($8)) t 'log(n_{terminals})', '' us 2:(log($9)) t 'log(n_{nonTerminals})'");
    //  gnuplot("load 'plt'");

  m /= double(K);
  v = ::sqrt(v/double(K) - m*m);
  v /= ::sqrt(double(K));
}

//===========================================================================

void ex1_HiddenTarget(){
  rai::String ex="hiddenTarget_";
  Problem_Options opt;
  CHECK_EQ(opt.problem, 2, "");
  ex <<opt.solver;
  if(opt.solver=="ELS") ex <<opt.wP <<opt.cP;
  ofstream fil(ex+".dat");
  for(uint dep=1;dep<=5;dep++){
    for(uint target=0;target<=5;target++){
      rai::setParameter<double>("Tree/maxDepth", dep);
      rai::setParameter<double>("Tree/hiddenTarget", target);
      double m, v;
      evalTimeToTarget(10, m, v);
      LOG(0) <<"depth: " <<dep <<" target: " <<target <<" mean: " <<m <<" meanSdv: " <<v;
      fil <<"depth: " <<dep <<" target: " <<target <<" mean: " <<m <<" meanSdv: " <<v <<endl;
    }
    fil <<endl;
  }
}

//===========================================================================

void ex4_lgpSolver(){
  rai::String problem = rai::getParameter<rai::String>("problem", STRING("none"));
  uint evalLimit = rai::getParameter<double>("LGP/evalLimit");
  Problem_Options opt;
  rai::NodeGlobal opt2;
  
  rai::String key = problem.getSubString(problem.find('/', true)+1, -1);
  key <<'_' <<opt.solver;
  if(opt.solver=="ELS") key <<opt2.level_wP <<opt2.level_cP;

  LOG(0) <<"============== run: " <<key;

  ofstream filStop(STRING("ex4_data/" <<key <<".STOP.dat"));

  uint K=10;
  for(uint k=0;k<K;k++){
    rai::LGP_Tool lgp(problem);

//    LOG(0) <<"LGP info:";
//    lgp.report(cout);

    rai::AStar astar(lgp.lgproot);
    if(opt.solver=="TreePolicy"){ NIY; }//astar.mode=astar.treePolicy;
    else if(opt.solver=="RoundRobin") astar.mode=astar.FIFO;

    rai::LGPComp_root* root = lgp.lgproot.get();
    ofstream fil(STRING("ex4_data/" <<key <<'.' <<k <<".dat"));
    uint c_tot=0;
    bool solFound=false;
    while(c_tot<evalLimit && astar.queue.N){
      astar.step();

      //check for solutions
      uint solutions=0;
      for(rai::TreeSearchNode *n:astar.solutions) if(n->isFeasible) solutions++;
      if(!solFound && solutions>0){ //write first time solution is found
	filStop <<astar.steps <<' ' <<c_tot <<' ' <<root->c_tot <<' '
		<<astar.mem.N <<' ' <<astar.solutions.N <<' ' <<solutions <<endl;
	solFound=true;
      }

      //check for compute increment
      if(floor(root->c_tot)>c_tot){ //write into the file only when a new c_tot integer is reached -- for comparability
        c_tot++;
        fil <<astar.steps <<' ' <<c_tot <<' ' <<root->c_tot <<' '
	    <<astar.mem.N <<' ' <<astar.solutions.N <<' ' <<solutions <<endl;
        cout <<"== steps: " <<astar.steps <<" c_tot: " <<c_tot <<' ' <<root->c_tot <<" nodes: " <<astar.mem.N <<" terminals: " <<astar.solutions.N <<" solutions: " <<solutions <<endl;
      }
    }

    if(!solFound){
      filStop <<astar.steps <<' ' <<c_tot <<' ' <<root->c_tot <<' '
	      <<astar.mem.N <<' ' <<astar.solutions.N <<' ' <<0 <<endl;
    }
  }

  procEval(STRING("ex4_data/" <<key <<'.'), K);
  rai::system(STRING("cp ex4_data/" <<key <<'.' <<"VAR.dat ex4_data/LGP_" <<key <<".dat"));
}

//===========================================================================

void ex4_all(){
  StringA problems = rai::getParameter<StringA>("problems");
  for(rai::String& p: problems){
    rai::setParameter<rai::String>("problem", p);
    ex4_lgpSolver();
  }
}

//===========================================================================

int main(int argn, char** argv){
  rai::initCmdLine(argn, argv);


  rnd.seed(1);
//  rnd.clockSeed();


//  ex1_HiddenTarget();
//  ex2_Random(); return 0;
//  ex3_bandit(); return 0;
  //ex4_lgpSolver(); return 0;
  ex4_all(); return 0;
  return 0;

  auto root = make_shared<HardBandit>(nullptr,0);
//  auto root = make_shared<HardBandit>(nullptr,0);

  rai::AStar astar(root);
  root->opt->fil = make_shared<ofstream>("z.dat");

  for(uint k=0;k<100;k++){
    astar.run(100);
    printTree(astar.mem);
    gnuplot("plot 'z.dat' us 2:3 t 'y' w p, '' us 2:4 t 'y_{mean}' lw 2, '' us 2:7 t 'y_{max}' lw 2, '' us 2:(log($8)) t 'log(n_{terminals})', '' us 2:(log($9)) t 'log(n_{nonTerminals})'");
    rai::wait();
  }

  return 0;
}
