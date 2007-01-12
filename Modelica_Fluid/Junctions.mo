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
    
    Modelica_Fluid.Interfaces.FluidPort_b port_1(redeclare package Medium = 
          Medium) annotation (extent=[-120,-10; -100,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
          Medium) annotation (extent=[100,-10; 120,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_3(redeclare package Medium = 
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
    Modelica_Fluid.Interfaces.FluidPort_b port_1(redeclare package Medium = 
          Medium) annotation (extent=[-120,-10; -100,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
          Medium) annotation (extent=[100,-10; 120,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_3(redeclare package Medium = 
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
  
end Junctions;
