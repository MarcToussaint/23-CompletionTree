#include "solver_uct.h"
#include <math.h>

void UCT_Solver::round(uint maxSteps){
    uint steps=0;
    simulator.reset();
    double totalReturn = 0;

    UCT_Node *node_v = root.get();
    uintA actions;
    Domain::StepReturn ret;

    //-- TREE POLICY
    while(steps < maxSteps                    //# we have compute budget
          && node_v->children.N == actions.N  //# we are 'inside' the tree: children for each action -> UCB to select the most promising
          && !ret.terminal){                  //# we're not at a terminal yet

        // compute the UCB scores for all children
        arr scores(node_v->children.N);
        uint i=0;
        for(auto& child: node_v->children){
            scores(i++) = child->data_Q / child->data_n + beta * sqrt(2. * ::log(node_v->data_n)/child->data_n);
        }
                
        // pick the child with highest
        node_v = node_v->children(argmax(scores)).get();

        // execute that action
        ret = simulator.step(node_v->action);
        totalReturn += ret.reward;
        steps++;
    }

    //-- ROLLOUT
    if(!ret.terminal){
        CHECK(node_v->children.N < actions.N, "");
        //# we have less children than actions -> we'll have to choose a new action, create a new child node, && preform a random rollout

        //# collect all fresh actions
        uintA virginActions = actions.copy();
        for(auto& child:node_v->children){
            virginActions.remove(child->action);
        }

        //# choose randomly from the fresh ones
        uint a = virginActions.rndElem();

        //# create a new child node
        auto newChild = make_shared<UCT_Node>(node_v, a);
        node_v->children.append(newChild);
        node_v = newChild.get();

        //# finish a rollout as long as we have budget
        while(!ret.terminal){
            ret = simulator.step(a);
            totalReturn += ret.reward;
            steps++;
            a = actions.rndElem();
        }
    }

    //-- BACKUP
    if(ret.terminal){ // (in case we couldn't finish the rollout)
        while(node_v){
            node_v->data_Q += totalReturn;
            node_v->data_n += 1;
            node_v = node_v->parent;
        }
    }
}

