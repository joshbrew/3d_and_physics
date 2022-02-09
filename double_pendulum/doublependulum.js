//Classic double pendulum complex system simulation using lagrangian mechanics

/* Basic complex system example...


The Double Pendulum:
___________
     |
     |
     |
 <<- O   ?
      \
       \
        \
     ?   O ->>

Variables:

g     gravity 9.81 m/s^2
m1    mass 1 (mid)
m2    mass 2 (end)
L1    length 1 (top->mid)
L2    length 2 (mid->end)
Θ1    angle 1 (top->mid)
Θ1'   angular velocity 1 
Θ1''  angular accel 1
Θ2    angle 2 (mid->end)
Θ2'   angular velocity 2
Θ2''  angular accel 2

Given any initial conditions for each variable (Θ1', Θ2', Θ1'', Θ2'' can be zero):

Θ1'' = ( -g*(2*m1 + m2)*sin(Θ1) - m2*g*sin(Θ1 - 2*Θ2) - 2*sin(Θ1 - Θ2)*m2*(Θ1'^2 * L2 + Θ1'^2 * L1 * cos(Θ1 - Θ2)) ) 
            /   L1*(2*m1 + m2 - m2*cos(2*Θ1 - 2*Θ2))

Θ2'' = ( 2*sin(Θ1-Θ2)*(Θ1'^2 * L1*(m1 + m2) + g*(m1 + m2)*cos(Θ1) + Θ2'^2 * L2*m2*cos(Θ1 - Θ2)) )
            /   L2*(2*m1 + m2 - m2*cos(2*Θ1 - 2*Θ2))

Where:
Θ1' = (Θ1(t_n) - Θ1(t_n-1))  /  (t_n - t_n-1) = Θ1''*dt
Θ2' = (Θ2(t_n) - Θ2(t_n-1))  /  (t_n - t_n-1) = Θ2''*dt

So calc first step with velocities and accelerations at whatever initial values.
Then just step forward from there by t.

After creating a path, let's see if we can predict it using the euler-lagrange equations and then 
minimize for the 'correct' initial constraints via maxent using some minor variations on initial conditions


*/

//Double Pendulum class
//le platonic ideal
export class DPend {
    
    g=9.81 //m/s^2
    m1=1 //mass 1
    m2=1 //mass 2
    L1=1
    L2=1
    Θ1=DPend.dtr(30)
    Θ2=DPend.dtr(30)
    Θ1dt=0
    Θ1ddt=0
    Θ2dt=0
    Θ2ddt=0
    x1=0 //set on first step by initial angles, lengths etc.
    y1=0
    x2=0
    y2=0
    loss=0.99

    //can set the initial values
    constructor(initialConditions={
        g:9.81, //m/s^2
        m1:1, //mass 1
        m2:1, //mass 2
        L1:1,
        L2:1,
        Θ1:DPend.dtr(-30),
        Θ2:DPend.dtr(30),
        Θ1dt:0,
        Θ1ddt:0,
        Θ2dt:0,
        Θ2ddt:0,
        loss:0.99 //loss multiplier
    }) {
        Object.assign(this,initialConditions);
        this.step(0); //set x and y coordinates from initial conditions
    }

    //degree to radian
    static dtr(d) {
        return Math.PI*d/180;
    }

    static rtd(r) {
        return r*180/Math.PI;
    }

    step(dt) {

        const t = this;

        //update angular accelerations
        t.Θ1ddt = ( -t.g*(2*t.m1 + t.m2)*Math.sin(t.Θ1) - t.m2*t.g*Math.sin(t.Θ1 - 2*t.Θ2) - 2*Math.sin(t.Θ1 - t.Θ2)*t.m2*(t.Θ1dt*t.Θ1dt * t.L2 + t.Θ1dt*t.Θ1dt * t.L1 * Math.cos(t.Θ1 - t.Θ2)) ) /   
                        (t.L1*(2*t.m1 + t.m2 - t.m2*Math.cos(2*t.Θ1 - 2*t.Θ2)) );

        t.Θ2ddt = ( 2*Math.sin(t.Θ1-t.Θ2)*(t.Θ1dt*t.Θ1dt * t.L1*(t.m1 + t.m2) + t.g*(t.m1 + t.m2)*Math.cos(t.Θ1) + t.Θ2dt*t.Θ2dt * t.L2*t.m2*Math.cos(t.Θ1 - t.Θ2)) ) /   
                        (t.L2*(2*t.m1 + t.m2 - t.m2*Math.cos(2*t.Θ1 - 2*t.Θ2)) );
    
        if(dt !== 0) {
            //update angular velocities
            t.Θ1dt = t.Θ1dt*t.loss + t.Θ1ddt*dt// - dt*t.g*Math.abs(Math.sin(t.Θ1));
            t.Θ2dt = t.Θ2dt*t.loss + t.Θ2ddt*dt// - dt*t.g*Math.abs(Math.sin(t.Θ2));

            //update angles
            t.Θ1=t.Θ1+t.Θ1dt*dt;
            t.Θ2=t.Θ2+t.Θ2dt*dt;
        }
        //update coordinates
        t.x1 = t.L1*Math.sin(t.Θ1);
        t.x2 = t.x1 + t.L2*Math.sin(t.Θ2);

        t.y1 = -t.L1*Math.cos(t.Θ1);
        t.y2 = t.y1 - t.L2*Math.cos(t.Θ2);

        return true;
    }
                
    //https://diego.assencio.com/?index=1500c66ae7ab27bb0106467c68feebc6
    lagrangian() {
        const t = this;
        
        let KineticEn = 0.5*t.m1*t.L1*t.L1 * t.Θ1dt*t.Θ1dt
                            + 0.5*t.m2*(t.L1*t.L1 * t.Θ1dt*t.Θ1dt
                            + 2*t.L1*t.L2*t.Θ1dt*t.Θ2dt*Math.cos(t.Θ1-t.Θ2));
        let PotentialEn = -(t.m1 + t.m2)*t.g*t.L1*Math.cos(t.Θ1) 
                            - t.m2*t.g*t.L2*Math.cos(t.Θ2);
        
        let L = KineticEn - PotentialEn;

        return L;
    }

    hamiltonian() {
        const t = this;
        
        let KineticEn = 0.5*t.m1*t.L1*t.L1 * t.Θ1dt*t.Θ1dt 
                            + 0.5*t.m2*(t.L1*t.L1 * t.Θ1dt*t.Θ1dt
                            + 2*t.L1*t.L2*t.Θ1dt*t.Θ2dt*Math.cos(t.Θ1-t.Θ2));
        let PotentialEn = -(t.m1 + t.m2)*t.g*t.L1*Math.cos(t.Θ1) 
                            - t.m2*t.g*t.L2*Math.cos(t.Θ2);
        
        let H = KineticEn + PotentialEn;

        return H;
    }

    canonical_momenta() {
        const t = this;

        let ðLðΘ1 = (t.m1+t.m2)*t.Θ1dt*t.L1*t.L1
                        + t.m2*t.L1*t.L2*t.Θ2dt*Math.cos(t.Θ1 - t.Θ2);

        let ðLðΘ2 = t.m2*t.Θ2dt*t.L2*t.L2
                        + t.m2*t.L1*t.L2*t.Θ1dt*Math.cos(t.Θ1 - t.Θ2);

        return {
            p1:ðLðΘ1,  //p = m*v ... so v = p/m and x = (p/m)*dt
            p2:ðLðΘ2
        };
    }

    //step forward using the canonical momenta obtained from the euler-lagrangian equations
    lstep(dt) {
        const t = this; 

        let p = this.canonical_momenta();

        t.Θ1dt = t.Θ1dt*t.loss + dt*p.p1/t.m1 - dt*t.g*Math.abs(Math.sin(t.Θ1));
        t.Θ2dt = t.Θ2dt*t.loss + dt*p.p2/t.m2 - dt*t.g*Math.abs(Math.sin(t.Θ2));

        t.Θ1 += t.Θ1dt*dt; 
        t.Θ2 += t.Θ2dt*dt;
        
        //update coordinates
        t.x1 = t.L1*Math.sin(t.Θ1);
        t.x2 = t.x1 + t.L2*Math.sin(t.Θ2);

        t.y1 = -t.L1*Math.cos(t.Θ1);
        t.y2 = t.y1 - t.L2*Math.cos(t.Θ2);
    }

    alphas() {
        const t = this;

        let a1 = (t.L2/t.L1)*(t.m2/(t.m1+t.m2))*Math.cos(t.Θ1-t.Θ2);
        let a2 = (t.L1/t.L2)*Math.cos(t.Θ1-t.Θ2);

        return {a1, a2};

    }

    second_order_diffs() {
        const t = this;

        let f1 = -(t.L2/t.L1)*(t.m2/(t.m1+t.m2))*Math.sin(t.Θ1-t.Θ2)*t.Θ1dt*t.Θ1dt - (t.g/t.L1)*Math.sin(t.Θ1);
        let f2 = (t.L1/t.L2)*Math.sin(t.Θ1-t.Θ2) - (t.g/t.L2)*Math.sin(t.Θ2);

        return {f1, f2};
    }

    det() {
        let alphas = this.alphas();

        let det = 1-(alphas.a1*alphas.a2);

        return det;
    }

    //first order diff eq solution based on angle and angular velocity
    fstep(dt) {
        const t = this;
        
        let f = this.second_order_diffs();
        let a = this.alphas();

        t.Θ1dt = t.Θ1dt*t.loss - dt*t.g*Math.abs(Math.sin(t.Θ1)) + dt*(f.f1 - a.a1*f.f2) / (1 - a.a1*a.a2);
        t.Θ2dt = t.Θ2dt*t.loss - dt*t.g*Math.abs(Math.sin(t.Θ2)) + dt*(f.f2 - a.a1*f.f1) / (1 - a.a1*a.a2);

        t.Θ1 += t.Θ1dt*dt; 
        t.Θ2 += t.Θ2dt*dt;

        //update coordinates
        t.x1 = t.L1*Math.sin(t.Θ1);
        t.x2 = t.x1 + t.L2*Math.sin(t.Θ2);

        t.y1 = -t.L1*Math.cos(t.Θ1);
        t.y2 = t.y1 - t.L2*Math.cos(t.Θ2);
        
    }

    //given the same initial conditions, all 3 methods should perform the same


}



//now draw the path of each ball in the pendulum using the lagrangian




/*

Lagrangian mechanics:

L = T - V

Lagrangian = Kinetic Energy - Potential Energy (Change in energy)

e.g. harmonic oscillator:
L = mass * velocity^2 / 2 - spring_const * position^2 / 2

Alternatively: Hamiltonian = Kinetic Energy + Potential Energy (Total energy)

The path of the change in energy follows with the principle of least action:

S = integral( (T(q,dq/dt) - V(q))  dt  ) from t0 to t1.

Producing the constraining Euler-Lagrange equation:

ðL/ðq - d/dt * ðL/ðq = 0 

Or for three independent variables:

ðL/ðq - ð/ðx * ðL/ðq_x - ð/ðy * ðL/ðq_y - ð/ðz * ðL/ðq_z = 0 

If there are many possible constraints i.e. we are unsure of our model, 
the best constraints can often be calculated with the principle of maximum entropy:

H = - integral ( p ln(p/q) ) dt from t0 to t1 and the action is zero for each statistical moment

Where p are the observed values at time t and q are the predicted values.

This is because nature will want to maximize its energy dissipation i.e. thermodynamic entropy.
The same holds generally true in information systems which are still fundamentally energy-based.

*/



