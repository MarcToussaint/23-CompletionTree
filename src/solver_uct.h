#include <Logic/treeSearchDomain.h>
#include <Core/array.h>

#include <memory>
using std::shared_ptr;

struct Domain{
    struct StepReturn {
        double reward=0.;
        double observation=0.;
        bool terminal=false;
    };

    virtual void reset() = 0;
    virtual uint getNumDecisions() = 0;
    virtual StepReturn step(uint decision) = 0;
};

struct UCT_Node{
    UCT_Node *parent=0;
    int action=-1;
    rai::Array<shared_ptr<UCT_Node>> children;
    double data_Q=0.;
    uint data_n=0;

    UCT_Node(UCT_Node *_parent, int _action) : parent(_parent), action(_action) {}
};

struct UCT_Solver{
    Domain& simulator;
    double beta=1.;

    UCT_Solver(Domain& _P) : simulator(_P) {}
    shared_ptr<UCT_Node> root;
    void round(uint maxSteps=(1<<20));
};
