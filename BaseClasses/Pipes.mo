package Pipes 
partial model Flow1D_FV 
    import Modelica_Fluid.Types;
    import Modelica.Constants.*;
  replaceable package Medium = PackageMedium 
    extends Modelica.Media.Interfaces.PartialMedium "Fluid medium model" 
   annotation (choicesAllMatching=true);
    
//Discretization
  parameter Integer n(min=1) "Number of pipe segments";
  final parameter Integer np=if lumped_dp then 1 else n + 1 
      "Number of momentum balances"                                                     annotation(Dialog(tab="Advanced"),Evaluate=true);
    
//Advanced model options
  parameter Boolean allowFlowReversal=true 
      "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                     annotation(Dialog(tab="Advanced", group="Mass and energy balances"));
  parameter Boolean static=true "= true, no mass or energy is stored" 
                                annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not dynamicTerm),Evaluate=true);
  parameter Boolean lumped_dp=true 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                                                                              annotation(Dialog(tab="Advanced", group="Momentum balance", enable=not dynamicTerm),Evaluate=true);
  parameter Boolean kineticTerm=false " = true, include kinetic term" 
                                              annotation(Dialog(tab="Advanced", group="Momentum balance", enable=not lumped_dp),Evaluate=true);
  parameter Boolean dynamicTerm=false 
      " = true, include dynamic term, only if not lumped_dp and not static"                               annotation(Dialog(tab="Advanced", group="Momentum balance", enable=(not static and not lumped_dp)),Evaluate=true);
    
//Initialization
    parameter Types.InitWithGlobalDefault.Temp initOption=Types.
        InitWithGlobalDefault.UseGlobalFluidOption "Initialization option" 
    annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
      "Start value of pressure at port a" 
    annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
      "Start value of pressure at port b" 
    annotation(Dialog(tab = "Initialization"));
    final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
        p_a_start - (p_a_start - p_b_start)/(2*n),
        p_b_start + (p_a_start - p_b_start)/(2*n),
        n) else {(p_a_start + p_b_start)/2} "Start value of pressure";
    parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
              Medium.temperature_phX(
        (p_a_start + p_b_start)/2,
        h_start,
        X_start) "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
        Medium.specificEnthalpy_pTX(
        (p_a_start + p_b_start)/2,
        T_start,
        X_start) else Medium.h_default "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.MassFlowRate mflow_start "Start value for mass flow rate" 
                                                                             annotation(Evaluate=true, Dialog(tab = "Initialization"));
    final parameter SI.Pressure dp_start=p_a_start - p_b_start;
    
//Geometry parameters
  parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
      Modelica_Fluid.Types.CrossSectionTypes.Circular 
      "Type of cross section of pipe" 
    annotation (Dialog(tab="General", group="Pipe geometry"));
  parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Pipe geometry"));
  parameter SI.Diameter d_inner "Inner diameter of circular pipe" annotation(Dialog(group="Pipe geometry", enable=crossSectionType==1));
  parameter SI.Length h_inner "Inner height of rectangular pipe"    annotation(Dialog(group="Pipe geometry", enable=crossSectionType==2));
  parameter SI.Length w_inner "Inner width of rectangular pipe"    annotation(Dialog(group="Pipe geometry", enable=crossSectionType==2));
  parameter SI.Length P_inner=if crossSectionType == 1 then Modelica.Constants.pi*d_inner else if crossSectionType == 2 then 2*h_inner + 2*
      w_inner else 1 "Inner perimeter" 
                                      annotation(Dialog(tab="General", group="Pipe geometry", enable=crossSectionType==3));
  inner parameter SI.Area A_inner=if crossSectionType == 1 then Modelica.Constants.pi*d_inner*d_inner/4 else if crossSectionType
       == 2 then h_inner*w_inner else 1 "Inner cross section area" 
                                          annotation(Dialog(tab="General", group="Pipe geometry", enable=crossSectionType==3));
  final parameter SI.Volume V=A_inner*length "Volume" 
                                                     annotation(Dialog(tab="General", group="Pipe geometry"));
    
//Pressure Drop
  replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.QuadraticTurbulent 
    extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
      "Characteristic of wall friction" 
                                       annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
  parameter SI.Diameter d_h=4*A_inner/P_inner "Hydraulic diameter" annotation(Dialog(tab="General", group="Pressure loss"));
  parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Evaluate=true);
  parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
  parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
  parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                            annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                             annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
  parameter Boolean show_Re=false 
      "= true, if Reynolds number is included for plotting" 
     annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter Boolean from_dp=true 
      " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
    annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
  parameter SI.AbsolutePressure dp_small=1 
      "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=from_dp and WallFriction.use_dp_small));
  parameter SI.MassFlowRate m_flow_small=0.01 
      "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
    annotation(Dialog(tab="Advanced",group="Pressure loss", enable=not from_dp and WallFriction.use_m_flow_small));
  SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
      m_flow,
      (eta_a + eta_b)/2,
      diameter) if                                                                                               show_Re 
      "Reynolds number of pipe flow";
    
//Total quantities
  SI.Energy[n] U "Internal energy of fluid";
  SI.Mass[n] m "Mass of fluid";
  SI.Mass[n,Medium.nXi] mXi "Masses of independent components in the fluid";
    
//Flow quantities
  inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -inf else 
              0, each start=mflow_start, each fixed=false) 
      "Mass flow rates of fluid across segment boundaries";
  SI.Velocity[n+1] v "velocity at volume boundaries (for display purposes)";
  Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
      "Independent mass flow rates across segment boundaries";
  Medium.EnthalpyFlowRate[n + 1] H_flow 
      "Enthalpy flow rates of fluid across segment boundaries";
  Real[n+1] I_flow;
  SI.Pressure[np] dp(start=dp0) 
      "Pressure drop due to friction loss and gravity";
    
//Source terms
  Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
  Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
      "Independent mass flow rates, source or sink";
  SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    
//Fluid ports
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
        Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
      "Fluid inlet port" 
                       annotation (extent=[-112,-10; -92,10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
        Medium, m_flow(max=if allowFlowReversal and not static then +inf else 0)) 
      "Fluid outlet port" 
                        annotation (extent=[90,-10; 110,10]);
  Medium.BaseProperties[n] medium(
    each preferredMediumStates=if static then false else true,
    p(start=p_start),
    each h(start=h_start),
    each T(start=T_start),
    each Xi(start=X_start[1:Medium.nXi]));
 /* Medium.BaseProperties medium_a(
    preferredMediumStates=false,
    p(start=p_a_start),
    h(start=h_start),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));
  Medium.BaseProperties medium_b(
    preferredMediumStates=false,
    p(start=p_b_start),
    h(start=h_start),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));*/
  annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
          color=69,
          gradient=2,
          fillColor=69))),
      Documentation(info="<html>
<p>
From Katrins email, Nov. 28, 2005:
</p>
 
<p>
Distributed volume model, properties and flow variables are arrays, no components as in the isolated pipe. Momentum and energy balances on a staggered grid, half a momentum balance on each end of the pipe. The medium properties in the ports are those of the upstream volume. I am strongly in favour with not using the energy balance 2 (the one where the momentum balance has been substracted) here, because you are loosing all the benefits of a staggered grid. You need twice as many momentum balances to calculate the algebraic pressures at the volume boundary which appear now in the energy balance. (And I am not sure if this can be  properly handled by the tool). The pressure drop is then also part of the energy balance and needs to be in accordance with the chosen grid. I agree with you that neglecting potential and kinetic energy in this case might not comply with a highly accurate formulation for teaching purposes, but for most applications it is more than sufficient. However, velocity and gravity can play a significant role in the momentum balance, which should have the option to include those terms (-> dynamic pressure). Not intertwining the two balances has also the advantage to be able to neglect specific terms in one of the balances and not in both.
</p>
 
<p>
The model contains source terms in mass and energy balances, which are not determined here. Therefore it is a partial model and could also be used for reactions or partial condensing gases with neglectable liquid volume (-> i.e. moist air).
</p>
 
<pre>
Modelling options (via boolean flags) are:
- static or dynamic (mass and energy) balances
- lumped pressure drop
- lumped composition (not sure yet if that is feasible)
- including velocity and gravity term in momentum balance, perhaps also the dynamic term.
</pre>
 
<p>
One issue not solved yet: For pressure drop and velocity term the densities at the ports are required. A medium function computing density from p and h would be most convenient. 
</p>
</html>"));
    
  protected 
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[1]);//approximation
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[n]);//approximation
  SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else 
            (if use_eta_nominal then fill(eta_nominal, n) else 
      Medium.dynamicViscosity(medium));
  SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
  SI.Density d_a=if use_d_nominal then d_nominal else medium[1].d;//approximation
  SI.Density d_b=if use_d_nominal then d_nominal else medium[n].d;//approximation
  /*SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else noEvent(if (dp[1] >= 0) then 
      Medium.dynamicViscosity(medium_a) else Medium.dynamicViscosity(medium[1])));
  SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if 
      use_eta_nominal then eta_nominal else noEvent(if (dp[np] >= 0) then 
      Medium.dynamicViscosity(medium[n]) else Medium.dynamicViscosity(medium_b)));
  SI.Density d_a=if use_d_nominal then d_nominal else noEvent(if (dp[1] >= 0) then 
            medium_a.d else medium[1].d);
  SI.Density d_b=if use_d_nominal then d_nominal else noEvent(if (dp[np] >= 0) then 
           medium[n].d else medium_b.d);
 The outlet port medium model (medium_a or medium_b) describes the medium properties just after the outlet connection
  point which would give wrong densities at the pipe outlet for pressure drop calculations
  if fluid flows with substantially different densities are mixed*/
  outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
  parameter Types.Init.Temp initOption2=if initOption == Types.
      InitWithGlobalDefault.UseGlobalFluidOption then fluidOptions.
      default_initOption else initOption 
      annotation(Evaluate=true, Hide=true);
  parameter SI.Pressure[np] dp0={if lumped_dp then dp_start else dp_start/(if i
       > 1 and i < np then n else 2*n) for i in 1:np};
    
//Momentum balance terms
  //SI.Force[np] DI_flow "Delta momentum flow across flow grid boundaries";
  SI.Force[np] F_f "Friction force";
  SI.Force[np] F_p "Pressure forces";
    
initial equation 
  // Initial conditions
  if not static then
    if initOption2 == Types.Init.NoInit then
    // no initial equations
    elseif initOption2 == Types.Init.SteadyState then
    //steady state initialization
      if use_T_start then
        der(medium.T) = zeros(n);
      else
        der(medium.h) = zeros(n);
      end if;
      if not (lumped_dp or Medium.singleState) then
        der(medium.p) = zeros(n);
      elseif lumped_dp then
      //  der(medium[1].p) = 0;
      end if;
      if dynamicTerm then
        der(m_flow[1:n]) = zeros(n);
      end if;
      for i in 1:n loop
        der(medium[i].Xi) = zeros(Medium.nXi);
      end for;
    elseif initOption2 == Types.Init.InitialValues then
    //Initialization with initial values
      if use_T_start then
        medium.T = ones(n)*T_start;
      else
        medium.h = ones(n)*h_start;
      end if;
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
    //Steady state initialization for hydraulic states (p, m_flow)
      if use_T_start then
        medium.T = ones(n)*T_start;
      else
        medium.h = ones(n)*h_start;
      end if;
      if not (lumped_dp or Medium.singleState) then
        der(medium.p) = zeros(n);
      elseif lumped_dp then
        der(medium[1].p) = 0;
      end if;
      if dynamicTerm and not Medium.singleState then
        der(m_flow[1:n])=zeros(n);
      end if;
    else
      assert(false, "Unsupported initialization option");
    end if;
  end if;
    
equation 
    
  //Port medium models
/*  port_a.p = medium_a.p;
  port_b.p = medium_b.p;
  port_a.h = medium_a.h;
  port_b.h = medium_b.h;
  port_a.Xi = medium_a.Xi;
  port_b.Xi = medium_b.Xi;*/
    
  // Boundary conditions
  port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
  port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
  port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
  port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
  port_a.m_flow = m_flow[1];
  port_b.m_flow = -m_flow[n + 1];
    
  // Distributed flow quantities
  for i in 2:n loop
    H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
    mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
    v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/A_inner;
    if kineticTerm then
      I_flow[i] = Modelica_Fluid.Utilities.regSquare(m_flow[i], delta=max(0.0001,
      0.01*mflow_start))/A_inner/noEvent(if m_flow[i]>=0 then d[i-1] else d[i]);
    else
      I_flow[i] = 0;
    end if;
  end for;
  H_flow[1] = port_a.H_flow;
  H_flow[n + 1] = -port_b.H_flow;
  mXi_flow[1, :] = port_a.mXi_flow;
  mXi_flow[n + 1, :] = -port_b.mXi_flow;
  v[1] = m_flow[1]/d_a/A_inner;
  v[n + 1] = m_flow[n + 1]/d_b/A_inner;
  if kineticTerm then
    I_flow[1] = Modelica_Fluid.Utilities.regSquare(m_flow[1], delta=max(0.0001,
      0.01*mflow_start))/A_inner/d_a;
    I_flow[n + 1] = Modelica_Fluid.Utilities.regSquare(m_flow[n + 1], max(
      0.0001, 0.01*mflow_start))/d_b/A_inner;
  else
    I_flow[1] = 0;
    I_flow[n + 1] = 0;
  end if;
    
  // Total quantities
  for i in 1:n loop
    m[i] =V/n*medium[i].d;
    mXi[i, :] = m[i]*medium[i].Xi;
    U[i] = m[i]*medium[i].u;
  end for;
    
  //Mass and energy balance
  for i in 1:n loop
    if static then
      0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
      zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
    else
      der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
      der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
    end if;
  end for;
  for i in 1:n loop
  assert(allowFlowReversal or (m_flow[i]>=0),"Flow reversal not allowed");
  end for;
    
//Pressure drop and gravity
 if from_dp and not WallFriction.dp_is_zero then
    if lumped_dp then
      m_flow[1] = WallFriction.massFlowRate_dp(dp[1] - height_ab*fluidOptions.g*(
        d_a + d_b)/2, d_a, d_b, eta_a, eta_b, length, d_h, roughness,
        dp_small);
    else
      m_flow[1] = WallFriction.massFlowRate_dp(dp[1] - height_ab/2/n*
        fluidOptions.g*(d_a + d[1])/2, d_a, d[1], eta_a, eta[1],
        length/n/2, d_h, roughness, dp_small);
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp(dp[i] - height_ab/n*
          fluidOptions.g*(d[i-1] + d[i])/2, d[i-1],
          d[i], eta[i - 1], eta[i], length/n, d_h, roughness,
          dp_small);
      end for;
      m_flow[n + 1] = WallFriction.massFlowRate_dp(dp[np] - height_ab/n/2*
        fluidOptions.g*(d[n] + d_b)/2, d[n], d_b, eta[n], eta_b,
        length/n/2, d_h, roughness, dp_small);
    end if;
  else
    if lumped_dp then
      dp[1] = WallFriction.pressureLoss_m_flow(m_flow[1], d_a, d_b, eta_a,
        eta_b, length, d_h, roughness, m_flow_small) + height_ab*
        fluidOptions.g*(d_a + d_b)/2;
    else
      dp[1] = WallFriction.pressureLoss_m_flow(m_flow[1], d_a, d[1], eta_a,
        eta[1], length/n/2, d_h, roughness, m_flow_small) + height_ab/2/n*
        fluidOptions.g*(d_a + d[1])/2;
      for i in 2:n loop
        dp[i] = WallFriction.pressureLoss_m_flow(m_flow[i], d[i-1], d[i], eta[i-1],
          eta[i], length/n, d_h, roughness, m_flow_small) + height_ab/n*
          fluidOptions.g*(d[i-1] + d[i])/2;
      end for;
      dp[np] = WallFriction.pressureLoss_m_flow(m_flow[np], d[n], d_b, eta[n],
        eta_b, length/n/2, d_h, roughness, m_flow_small) + height_ab/n/2*
        fluidOptions.g*(d[n] + d_b)/2;
    end if;
  end if;
    
//Momentum Balance
if lumped_dp then
    F_p[1] = (port_a.p - port_b.p)*A_inner;
    F_f[1] = -dp[1]*A_inner;
    zeros(np) = F_p + F_f;
    medium.p = ones(n)*(port_a.p + port_b.p)/2;
  else
    F_p[1] = (port_a.p-medium[1].p)*A_inner;
    F_f[1] = -dp[1]*A_inner;
    (if dynamicTerm then der(m_flow[1])*length/n/2 else 0) = F_p[1] + F_f[1] + (I_flow[1]-I_flow[2])/2;
    for i in 2:n loop
      F_p[i] = (medium[i-1].p-medium[i].p)*A_inner;
      F_f[i] = -dp[i]*A_inner;
      (if dynamicTerm then der(m_flow[i])*length/n else 0) = F_p[i] + F_f[i] + (I_flow[i-1]-I_flow[i+1])/2;
    end for;
    F_p[np] = (medium[n].p-port_b.p)*A_inner;
    F_f[np] = -dp[np]*A_inner;
    (if dynamicTerm then der(m_flow[n + 1])*length/n/2 else 0) = F_p[np] + F_f[np] + (I_flow[n]-I_flow[n+1])/2;
  end if;
    
end Flow1D_FV;
  
model PortVolume 
    "Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
  import Modelica_Fluid.Types;
  extends Interfaces.PartialInitializationParameters;
    
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
  parameter SI.Volume V "Volume";
    
  Interfaces.FluidPort_a port(
    redeclare package Medium = Medium) "Fluid port" 
    annotation (extent=[-10, -10; 10, 10], rotation=0);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
    annotation (extent=[-10,90; 10,110]);
    
  Medium.BaseProperties medium(preferredMediumStates=true,
              p(start=p_start), T(start=T_start),
              h(start=h_start), Xi(start= X_start[1:Medium.nXi]));
  SI.Energy U "Internal energy of fluid";
  SI.Mass m "Mass of fluid";
  SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
  annotation (
   Icon(
      Ellipse(extent=[-100, 100; 100, -100], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=3,
          fillColor=68,
          rgbfillColor={170,213,255})),
      Text(extent=[-150,-100; 150,-150], string="%name")),
                         Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport. This splitting is only possible under certain assumptions.
</p>
</html>"),
    Diagram);
equation 
  // medium properties set by port values
    port.p = medium.p;
    port.h = medium.h;
    port.Xi = medium.Xi;
    thermalPort.T = medium.T;
    
  // Total quantities
     m    = V*medium.d "Total Mass";
     mXi = m*medium.Xi "Independent component masses";
     U    = m*medium.u "Internal energy";
    
  // Mass and energy balance
     der(m)    = port.m_flow "Total mass balance";
     der(mXi)  = port.mXi_flow "Independent component mass balance";
     der(U)    = port.H_flow + thermalPort.Q_flow "Energy balance";
    
initial equation 
  // Initial conditions
  if initOption2 == Types.Init.NoInit then
    // no initial equations
  elseif initOption2 == Types.Init.InitialValues then
    if not Medium.singleState then
       medium.p = p_start;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  elseif initOption2 == Types.Init.SteadyState then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    der(medium.h) = 0;
    der(medium.Xi) = zeros(Medium.nXi);
  elseif initOption2 == Types.Init.SteadyStateHydraulic then
    if not Medium.singleState then
       der(medium.p) = 0;
    end if;
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    medium.Xi = X_start[1:Medium.nXi];
  else
    assert(false, "Unsupported initialization option");
  end if;
end PortVolume;

package HeatTransfer 
  partial model PartialPipeHeatTransfer 
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
    parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
    SI.HeatFlowRate[n] Q_flow "Heat flow rates";
    parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
    parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
        "Thermal port" 
      annotation (extent=[-20,60; 20,80]);
    SI.Temperature[n] T;
  equation 
      
    annotation (Icon(Ellipse(extent=[-60,64; 60,-56], style(
            color=42,
            rgbcolor={127,0,0},
            gradient=3,
            fillColor=1,
            rgbfillColor={232,0,0})), Text(
          extent=[-38,26; 40,-14],
          style(
            color=42,
            rgbcolor={127,0,0},
            gradient=3,
            fillColor=1,
            rgbfillColor={232,0,0},
            fillPattern=7),
          string="%name")), Documentation(info="<html>
Base class for heat transfer models that can be used in model <b>Pipe</b>.
</html>"));
  end PartialPipeHeatTransfer;
    
  model PipeHT_constAlpha 
    extends PartialPipeHeatTransfer;
    parameter SI.CoefficientOfHeatTransfer alpha0=200;
    annotation(structurallyIncomplete=true, Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient
</html>"));
  equation 
    for i in 1:n loop
      thermalPort[i].Q_flow=alpha0*A_h/n*(thermalPort[i].T-T[i]);
    end for;
    thermalPort.Q_flow=Q_flow;
  end PipeHT_constAlpha;
  annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));
end HeatTransfer;
end Pipes;
