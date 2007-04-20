package Junctions "Junction components" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  
  model Splitter 
    "Splitting/joining component with static balances for an infinitesimal control volume" 
    
    annotation(Documentation(info="<html>
  This model is the simplest implementation for a splitting/joining component for
  three flows. Its use is not required. It just formulates the balance
  equations in the same way that the connect symmantics would formulate them anyways.
  The main advantage of using this component is, that the user does not get
  confused when looking at the specific enthalpy at each port which might be confusing
  when not using a splitting/joining component. The reason for the confusion is that one exmanins the mixing
  enthalpy of the infinitesimal control volume introduced with the connect statement when
  looking at the specific enthalpy in the connector which
  might not be equal to the specific enthalpy at the port in the \"real world\".</html>"));
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    
    Modelica.Fluid.Interfaces.FluidPort_b port_1(redeclare package Medium = 
          Medium) annotation (extent=[-120,-10; -100,10]);
    Modelica.Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
          Medium) annotation (extent=[100,-10; 120,10]);
    Modelica.Fluid.Interfaces.FluidPort_b port_3(redeclare package Medium = 
          Medium) annotation (extent=[-10,100; 10,120]);
    
    Medium.AbsolutePressure p "Pressure";
    Medium.SpecificEnthalpy hMix "Mixing enthalpy";
    Medium.MassFraction Xi[Medium.nXi] 
      "Independent mixture mass fractions m_i/m";
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    annotation (Icon(Polygon(points=[-100,60; -60,60; -60,100; 60,100; 60,60; 100,
              60; 100,-60; -100,-60; -100,60], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=9,
            rgbfillColor={175,175,175},
            fillPattern=1)), Polygon(points=[-100,40; -40,40; -40,100; 40,100; 40,
              40; 100,40; 100,-40; -100,-40; -100,40], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=2,
            fillColor=69,
            rgbfillColor={0,128,255}))));
  equation 
    port_1.m_flow + port_2.m_flow + port_3.m_flow = 0 "Mass balance";
    port_1.mXi_flow + port_2.mXi_flow + port_3.mXi_flow = zeros(Medium.nXi) 
      "Component mass balances";
    port_1.mC_flow + port_2.mC_flow + port_3.mC_flow = zeros(Medium.nC) 
      "Trace substance mass balances";
    
    port_1.H_flow = semiLinear(port_1.m_flow, port_1.h, hMix);
    port_2.H_flow = semiLinear(port_2.m_flow, port_2.h, hMix);
    port_3.H_flow = semiLinear(port_3.m_flow, port_3.h, hMix);
    
    port_1.mXi_flow = semiLinear(port_1.m_flow, port_1.Xi, Xi);
    port_2.mXi_flow = semiLinear(port_2.m_flow, port_2.Xi, Xi);
    port_3.mXi_flow = semiLinear(port_3.m_flow, port_3.Xi, Xi);
    
    port_1.mC_flow = semiLinear(port_1.m_flow, port_1.C, C);
    port_2.mC_flow = semiLinear(port_2.m_flow, port_2.C, C);
    port_3.mC_flow = semiLinear(port_3.m_flow, port_3.C, C);
    
    // Momentum balances
    port_1.p = p;
    port_2.p = p;
    port_3.p = p;
    
    port_1.H_flow + port_2.H_flow + port_3.H_flow = 0 "Energy balance";
  end Splitter;
  
  annotation (Documentation(info="<html>
 
</html>"));
  model Junction_dynamic 
    "Splitting/joining component with static balances for a dynamic control volume" 
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Volume V=1e-5;
    SI.InternalEnergy U;
    SI.Mass m;
    SI.Mass[Medium.nXi] mXi;
    Modelica.Fluid.Interfaces.FluidPort_b port_1(redeclare package Medium = 
          Medium) annotation (extent=[-120,-10; -100,10]);
    Modelica.Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
          Medium) annotation (extent=[100,-10; 120,10]);
    Modelica.Fluid.Interfaces.FluidPort_b port_3(redeclare package Medium = 
          Medium) annotation (extent=[-10,100; 10,120]);
    
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium;
    annotation (Icon(Polygon(points=[-100,60; -60,60; -60,100; 60,100; 60,60; 100,
              60; 100,-60; -100,-60; -100,60], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=9,
            rgbfillColor={175,175,175},
            fillPattern=1)), Polygon(points=[-100,40; -40,40; -40,100; 40,100; 40,
              40; 100,40; 100,-40; -100,-40; -100,40], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=2,
            fillColor=69,
            rgbfillColor={0,128,255}))));
  equation 
    port_1.m_flow + port_2.m_flow + port_3.m_flow = der(m) "Mass balance";
    port_1.mXi_flow + port_2.mXi_flow + port_3.mXi_flow = der(mXi) 
      "Component mass balances";
    port_1.mC_flow + port_2.mC_flow + port_3.mC_flow = zeros(Medium.nC) 
      "Trace substance mass balances";
    
    port_1.H_flow = semiLinear(port_1.m_flow, port_1.h,medium.h);
    port_2.H_flow = semiLinear(port_2.m_flow, port_2.h, medium.h);
    port_3.H_flow = semiLinear(port_3.m_flow, port_3.h, medium.h);
    
    port_1.mXi_flow = semiLinear(port_1.m_flow, port_1.Xi, medium.Xi);
    port_2.mXi_flow = semiLinear(port_2.m_flow, port_2.Xi, medium.Xi);
    port_3.mXi_flow = semiLinear(port_3.m_flow, port_3.Xi, medium.Xi);
    
    port_1.mC_flow = semiLinear(port_1.m_flow, port_1.C, C);
    port_2.mC_flow = semiLinear(port_2.m_flow, port_2.C, C);
    port_3.mC_flow = semiLinear(port_3.m_flow, port_3.C, C);
    
    // Momentum balance
    //Suitable for compressible media
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;
    
    U=m*medium.u;
    mXi=m*medium.Xi;
    m=medium.d*V;
    
    port_1.H_flow + port_2.H_flow + port_3.H_flow = der(U) "Energy balance";
  end Junction_dynamic;
  
  model MassFlowRatio "simple flow multiplier" 
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
    parameter Real ratio=1 "flow multiplier from port a to port b";
    Modelica.Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium=Medium) 
      annotation (extent=[-110,-10; -90,10]);
    Modelica.Fluid.Interfaces.FluidPort_b port_b(
                                  redeclare package Medium = Medium) 
      annotation (extent=[90,-10; 110,10]);
  equation 
    port_b.m_flow = ratio*port_a.m_flow;
    port_b.H_flow = ratio*port_a.H_flow;
    port_b.mXi_flow = ratio*port_b.mXi_flow;
    port_b.mC_flow = ratio*port_b.mC_flow;
    
    port_a.p = port_b.p;
    port_a.H_flow = semiLinear(
      port_a.m_flow,
      port_a.h,
      port_b.h);
    port_a.mXi_flow = semiLinear(
      port_a.m_flow,
      port_a.Xi,
      port_b.Xi);
    port_a.mC_flow = semiLinear(
      port_a.m_flow,
      port_a.C,
      port_b.C);
    annotation (Icon(
        Line(points=[-80,0; 80,80], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,0], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,-80], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,40], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,-40], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4))), Documentation(info="<html>
This model reflects a simple flow multiplication, which is very helpful in cases where the flow is evenly distributed to several parallel flow paths which are identical in their dimensions and boundary conditions, as e.g. in heat exchangers. Only one of the parallel pipes needs to be simulated then. All flow variables in <b>port_b</b> are achieved by multiplying those in <b>port_a</b> with parameter <b>ratio</b>.
</html>"));
  end MassFlowRatio;
  
  model HeatFlowRatio "simple heat flow multiplier" 
    parameter Real ratio=1 "heat flow ratio from port_a to port_b";
    annotation (Icon(
        Line(points=[-80,0; 80,80], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,40], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,0], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,-40], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,-80], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4))), Documentation(info="<html>
Simple model for heat flow multiplication between the two ports. The heat flow rate in port_a is multiplied by parameter <b>ratio</b> to achieve the rate at port_b. Both temperatures are equal, the model may be used e.g. for parallel pipes in heat exchangers.
</html>"));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a 
      annotation (extent=[-100,-10; -80,10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b 
      annotation (extent=[80,-10; 100,10]);
  equation 
    port_b.Q_flow=ratio*port_a.Q_flow;
    port_b.T=port_a.T;
  end HeatFlowRatio;
end Junctions;
