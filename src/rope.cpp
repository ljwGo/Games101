#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

float damping_factor = 0.00005;

// reference: http://games-cn.org/forums/topic/guanyuzuoye8deyixiewentijieda/
namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        Vector2D offset = (end - start) / (num_nodes - 1);

        Mass* startMass = new Mass(start, node_mass, true);
        // startMass->velocity = Vector2D(0, 0);  default none argument construct x and y are zero
        // startMass->forces = Vector2D(0, 0);
        masses.push_back(startMass);

        Mass* previouMass = startMass;
        Mass* nextMass = startMass;
        Vector2D nextPos = start;

        for (int i=0; i<num_nodes-1; i++){
            // generate mass
            previouMass = nextMass;
            nextPos += offset;

            nextMass = new Mass(nextPos, node_mass, false);
            masses.push_back(nextMass);

            // generate spring
            Spring* newSpring = new Spring(previouMass, nextMass, k);
            springs.push_back(newSpring);
        }

//        Comment-in this part when you implement the constructor
//        for (auto &i : pinned_nodes) {
//            masses[i]->pinned = true;
//        }
    }

    // An error will cause to a trouble. inverse the spring. 
    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Mass* m1 = s->m1;
            Mass* m2 = s->m2;
            
            float springLength = (m1->position - m2->position).norm();
            float forceSize = s->k * -(springLength - s->rest_length);

            Vector2D m2Tom1 = m1->position - m2->position;
            Vector2D m1Tom2 = m2->position - m1->position;

            // pay attention to force direction between stretch and compress
            Vector2D m2Force = m1Tom2 / m1Tom2.norm() * forceSize;
            Vector2D m1Force = m2Tom1 / m2Tom1.norm() * forceSize;

            m1->forces += m1Force;
            m2->forces += m2Force;
        }

        for (auto &m : masses)
        {
            // pinned mass's position has be fixed.
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += gravity;

                // explicit euler method
                // m->position += m->velocity * delta_t;
                // Vector2D accelerate = m->forces / m->mass;
                // m->velocity += accelerate * delta_t;

                // halt-implicit euler method
                // Vector2D accelerate = m->forces / m->mass;
                // m->velocity += accelerate * delta_t;
                // m->position += m->velocity * delta_t;

                // TODO (Part 2): Add global damping
                // explicit euler damping force = -k_d * v;
                float const k_d = 0.01;
                Vector2D f_damping = m->velocity * -k_d;
                m->forces += f_damping;

                Vector2D accelerate = m->forces / m->mass;
                m->velocity += accelerate * delta_t;
                m->position += m->velocity * delta_t;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerletForce(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {   // constraints meaing is similar with checked in database.
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Mass* m1 = s->m1;
            Mass* m2 = s->m2;
            
            float springLength = (m1->position - m2->position).norm();
            float forceSize = s->k * -(springLength - s->rest_length);

            Vector2D m2Tom1 = m1->position - m2->position;
            Vector2D m1Tom2 = m2->position - m1->position;

            // pay attention to force direction between stretch and compress
            Vector2D m2Force = m1Tom2 / m1Tom2.norm() * forceSize;
            Vector2D m1Force = m2Tom1 / m2Tom1.norm() * forceSize;

            m1->forces += m1Force;
            m2->forces += m2Force;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                // explicit verlet: x(t+1) = x(t) + [x(t)-x(t-1)] + a * t * t;
                m->forces += gravity;
                Vector2D accelerate = m->forces / m->mass;

                // have no damping
                // m->position = m->position + m->position - m->last_position + accelerate * delta_t * delta_t;
                // m->last_position = temp_position;

                // TODO (Part 4): Add global Verlet damping
                m->position += (1 - damping_factor) * (m->position - m->last_position) + accelerate * delta_t * delta_t;
                m->last_position = temp_position;
            }
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Mass* m1 = s->m1;
            Mass* m2 = s->m2;
            
            float springLength = (m1->position - m2->position).norm();
            float displace = springLength - s->rest_length;

            Vector2D m2Tom1 = m1->position - m2->position;
            Vector2D m1Tom2 = m2->position - m1->position;

            // solution refer to website above.
            Vector2D m1Modifier = 0.5 * m1Tom2 / m1Tom2.norm() * displace;
            Vector2D m2Modifier = 0.5 * m2Tom1 / m2Tom1.norm() * displace;

            if (!m1->pinned){
                m1->position += m1Modifier;
            }
            
            if (!m2->pinned){
                m2->position += m2Modifier;
            }
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                // explicit verlet: x(t+1) = x(t) + [x(t)-x(t-1)] + a * t * t;
                Vector2D accelerate = gravity / m->mass;

                // have no damping
                // m->position = m->position + m->position - m->last_position + accelerate * delta_t * delta_t;
                // m->last_position = temp_position;

                // TODO (Part 4): Add global Verlet damping
                m->position += (1 - damping_factor) * (m->position - m->last_position) + accelerate * delta_t * delta_t;
                m->last_position = temp_position;
            }
        }
    }
}
