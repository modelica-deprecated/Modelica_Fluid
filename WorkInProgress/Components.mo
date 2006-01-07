package Components 
  
model IsolatedPipe 
    "Model of an isolated pipe consisting of n pipe segments/FiniteVolumes" 
    import SI = Modelica.SIunits;
    
  replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching = true);
    
  extends Modelica_Fluid.Interfaces.PartialInitializationParameters;
    
  parameter Integer nVolumes(min=1)=1 "Number of pipe segments/finite volumes";
    
  parameter SI.Length L "Length of pipe";
  parameter SI.AbsolutePressure dp_nominal(min=1.e-10) = 1 
      "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
    
  parameter SI.MassFlowRate m_flow_nominal = 1E-3 
      "Nominal mass flow rate at nominal pressure drop";
    
  parameter SI.Area A_a;
  parameter SI.Area A_b=A_a;
    
  parameter SI.Length Z_a=0;
  parameter SI.Length Z_b=Z_a;
    
  parameter Boolean dynamicMomentumBalance=false 
      "If false, der(m_flow) is neglected in momentum balance" 
                                                 annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeKineticTerm=false 
      "If false, d*v^2 is neglected in momentum balance" 
                                             annotation(Evaluate=true,
      Dialog(tab="Level of Detail"));
  parameter Boolean includeViscosity=false 
      "If false, artifical viscosity is neglected" 
                                          annotation(Evaluate=true, Dialog(tab=
          "Level of Detail"));
  parameter Real viscosityFactor1=0 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
  parameter Real viscosityFactor2=1 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
    
  Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare model Medium = Medium) 
              annotation (extent=[-120, -10; -100, 10]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare model Medium = Medium) 
              annotation (extent=[120, -10; 100, 10]);
  Modelica_Fluid.WorkInProgress.Utilities.PipeSegment pipeSegment[nVolumes](
      redeclare package Medium = Medium,
      each initOption = initOption,
      each p_start = p_start,
      each use_T_start = use_T_start,
      each T_start = T_start,
      each h_start = h_start,
      each X_start = X_start,
      each L=L/nVolumes,
      each dp_nominal=dp_nominal/nVolumes,
      each A_a=A_a "has to be corrected: linear distribution of A",
      each Z_a=Z_a "has to be corrected: linear distribution of Z",
      each m_flow_nominal=m_flow_nominal,
      each dynamicMomentumBalance=dynamicMomentumBalance,
      each includeKineticTerm=includeKineticTerm,
      each includeViscosity=includeViscosity,
      each viscosityFactor1=viscosityFactor1,
      each viscosityFactor2=viscosityFactor2);
    
annotation (Icon(
    Rectangle(extent=[-100, 60; 100, -60], style(color=0, fillColor=8)),
    Rectangle(extent=[-100, 34; 100, -36], style(
        color=69,
        gradient=2,
        fillColor=69)),
    Text(
      extent=[-150, 125; 150, 65],
      string="%name",
      style(gradient=2, fillColor=69)),
      Ellipse(extent=[-58,14; -28,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0})),
      Ellipse(extent=[22,14; 52,-14], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=0,
            rgbfillColor={0,0,0}))));
equation 
  connect(port_a, pipeSegment[1].port_a);
  connect(port_b, pipeSegment[nVolumes].port_b);
  for i in 1:nVolumes - 1 loop
    connect(pipeSegment[i].port_b, pipeSegment[i + 1].port_a);
  end for;
end IsolatedPipe;
  
model Pipe 
    
  extends Modelica_Fluid.WorkInProgress.Interfaces.Flow1D(
    Qs_flow=heat.Q_flow,
    ms_flow=zeros(n),
    msXi_flow=zeros(n, Medium.nXi),
    dp=friction.dp);
  parameter SI.Area area_h = P_inner*length "Heat transfer area" annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  parameter SI.Diameter d_h = 4*A_inner/P_inner "Hydraulic diameter" annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  parameter Boolean use_wall=false 
      "= true, use wall component between fluid and thermalPort" 
                                                                annotation(Dialog(tab="General", group="Wall - optional"),Evaluate=true);
  parameter SI.Diameter d_outer "Outer diameter of circular pipe" annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==1 and use_wall)));
  parameter SI.Length h_outer "Outer height of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length w_outer "Outer width of rectangular pipe"  annotation(Dialog(tab="General", group="Wall - optional", enable=(crossSectionType==2 and use_wall)));
  parameter SI.Length A_outer = if crossSectionType == 1 then Modelica.Constants.pi*d_outer*d_outer/4 else if crossSectionType == 2 then h_outer*w_outer else 1 
      "Outer cross section area" 
                               annotation(Dialog(tab="General", group="Wall - optional", enable=(use_wall and crossSectionType==3)));
  inner Medium.ThermodynamicState[n] state = medium.state;
    
  replaceable 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PartialPipeHeatTransfer(
    redeclare final package Medium = Medium,
    final n=n,
    final d_h=d_h,
    final A_h=area_h,
    T=medium.T) "Convective heat transfer" 
                annotation (Dialog(tab="Advanced", group="Heat transfer and friction loss"),choicesAllMatching, extent=[-20,-20;
        20,20]);
  replaceable 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_RoughSurface
      friction(
    redeclare final package Medium = Medium,
    final d_h = d_h,
    final length = length,
    final n=n,
    final np=np,
    final lumped_dp=lumped_dp,
    m_flow=if lumped_dp then ones(np)*m_flow[2] else m_flow,
    d_1 = if lumped_dp then ones(np)*d_port_a else {if i>1 then medium[i-1].d else d_port_a for i in 1:np},
    d_2 = if lumped_dp then ones(np)*d_port_b else {if i<np then medium[i].d else d_port_b for i in 1:np}) extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PartialPipeFriction(
    redeclare final package Medium = Medium,
    final d_h = d_h,
    final length = length,
    final n=n,
    final np=np,
    final lumped_dp=lumped_dp,
    m_flow=if lumped_dp then ones(np)*m_flow[2] else m_flow,
    d_1 = if lumped_dp then ones(np)*d_port_a else {if i>1 then medium[i-1].d else d_port_a for i in 1:np},
    d_2 = if lumped_dp then ones(np)*d_port_b else {if i<np then medium[i].d else d_port_b for i in 1:np}) 
      "Pressure drop due to friction" 
                                    annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"), choicesAllMatching, extent=[-58,-20;
        -18,20]);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
    annotation (extent=[-20,60; 20,80]);
  replaceable model Wall = 
        Modelica_Fluid.WorkInProgress.Components.Wall_constProps                    extends 
      Modelica_Fluid.WorkInProgress.Interfaces.PartialPipeWall "Wall model" 
                                                                   annotation(choicesAllMatching, Dialog(enable=use_wall, tab="General", group="Wall - optional"));
  Wall wall(final n=n, final a_inner=A_inner, final a_outer=A_outer, final 
        length=length, initOption=initOption) if use_wall 
                           annotation (extent=[10,20; 50,60]);
  annotation (Icon(Rectangle(extent=[-100,60; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)), Rectangle(extent=[-100,-40; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))),
                            Diagram,
      Documentation(info="<html>
<p>
From Katrins email, Nov. 28, 2005:
</p>
 
<p>
extends Interfaces.1DFlow. Pressure drop and heat transfer are added in terms of replaceable components. The main problem here is to make all required variables and parameters available to the respective component (medium state, pipe geometry, Medium functions, empirical parameters). Only those shared by all future replaceable models (the simple one parameter model and the highly sophisticated (fictitious) two phase Nusselt correlation) can be set by modifiers (which is not straightforward in Dymola at the moment if a contsraining clause is used).  Those not required by all models as i.e. viscosity and conductivitiy must be computed inside the component from medium properties made available via inner and outer. I always try to avoid this as it it as bit like free climbing, but in this case I see no better solution.
</p>
 
<p>
Martin, I have not tested your latest pressure drop implementation with this model, but will do so as soon as possible. However, it is used in a completely different way, that means as an array of components, not as a  base class, in order to be able to handle distributed flow. I will check if another implementation would be more practical.
</p>
 
<p>
The pipe model contains a Boolean flag useWall which determines if a wall component is added. Unfortunately the icon does not represent the difference. In this way a heat exchanger can be created using two instances of the pipe model, one with a wall and one without. If interested in transients it could also make sense to include a wall in an insulated pipe. 
</p>
 
</html>"));
equation 
if use_wall then
  connect(wall.thermalPort_a, thermalPort) annotation (points=[30,50; 30,60; 0,
          60; 0,70],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
  connect(wall.thermalPort_b, heat.thermalPort) 
                                              annotation (points=[30,30; 30,22;
          0,22; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
else
  connect(thermalPort, heat.thermalPort) 
                                       annotation (points=[0,70; 0,14],
      style(
      color=10,
      rgbcolor={95,95,95},
      fillColor=0,
      rgbfillColor={0,0,0},
      fillPattern=7));
end if;
end Pipe;
  
model HeatExchanger "Double pipe heat exchanger with neglectible outer wall" 
  //General
  parameter Integer n(min=1) "Spatial segmentation";
  replaceable package Medium_1 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Inner pipe medium" 
                                                    annotation(choicesAllMatching);
  replaceable package Medium_2 = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialMedium "Outer pipe medium" 
                                                    annotation(choicesAllMatching);
  parameter SI.Length di_1(min=0) "Inner diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length da_1(min=0) "Outer diameter of inner pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length di_2(min=0) "Inner diameter of outer pipe"     annotation(Dialog(tab="General", group="Dimensions"));
  parameter SI.Length length(min=0) "Length of both pipes" annotation(Dialog(tab="General", group="Dimensions"));
  //Wall
  parameter SI.Density d_wall "Density of wall material" annotation(Dialog(tab="General", group="Constant material properties"));
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material"                                      annotation(Dialog(tab="General", group="Constant material properties"));
  parameter SI.Temperature T_start_w "Start value of wall termperature" annotation(Dialog(tab="Initialization", group="Wall"));
  final parameter SI.Mass m_wall=sum(pipe_1.wall.m) "Wall mass";
  //Initialization pipe 1
  parameter Modelica_Fluid.Types.Init.Temp initOption_1 "Initialization option"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Boolean use_T_start_1=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.AbsolutePressure p_start_1=Medium_1.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Inner pipe"));
  parameter Medium_1.Temperature T_start_1=if use_T_start_1 then Medium_1.T_default else 
      Medium_1.temperature_phX(p_start_1, h_start_1, X_start_1) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = use_T_start_1));
  parameter Medium_1.SpecificEnthalpy h_start_1=if use_T_start_1 then 
      Medium_1.specificEnthalpy_pTX(p_start_1, T_start_1, X_start_1[1:Medium_1.nXi]) else Medium_1.h_default 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe", enable = not use_T_start_1));
  parameter Medium_1.MassFraction X_start_1[Medium_1.nX]=Medium_1.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Inner pipe", enable=(Medium_1.nXi > 0)));
  parameter Medium_1.MassFlowRate mflow_start_1 "Start value of mass flow rate"
                                    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Inner pipe"));
  //Initialization pipe 2
  parameter Modelica_Fluid.Types.Init.Temp initOption_2 "Initialization option"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Boolean use_T_start_2=true "Use T_start if true, otherwise h_start"
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.AbsolutePressure p_start_2=Medium_2.p_default 
      "Start value of pressure" 
    annotation(Dialog(tab = "Initialization", group = "Outer pipe"));
  parameter Medium_2.Temperature T_start_2=if use_T_start_2 then Medium_2.T_default else 
      Medium_2.temperature_phX(p_start_2, h_start_2, X_start_2) 
      "Start value of temperature" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = use_T_start_2));
  parameter Medium_2.SpecificEnthalpy h_start_2=if use_T_start_2 then 
      Medium_2.specificEnthalpy_pTX(p_start_2, T_start_2, X_start_2[1:Medium_2.nXi]) else Medium_2.h_default 
      "Start value of specific enthalpy" 
    annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe", enable = not use_T_start_2));
  parameter Medium_2.MassFraction X_start_2[Medium_2.nX]=Medium_2.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", group = "Outer pipe", enable=Medium_2.nXi>0));
  parameter Medium_2.MassFlowRate mflow_start_2 "Start value of mass flow rate"
                                       annotation(Evaluate=true, Dialog(tab = "Initialization", group = "Outer pipe"));
  //Advanced
  parameter Boolean lumped_dp 
      " = true, lumped pressure drop, reduces number of pressure states to one"
                              annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean static "= true, use quasistatic mass and energy balances" 
                           annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Boolean kineticTerm 
      " = true, include kinetic term in momentum balance" 
                                annotation(Evaluate=true, Dialog(tab="Advanced", group="Conservation of mass, energy, momentum"));
  parameter Real K1=1 "Enhancement factor for outer pipe heat transfer area" 
                                                                          annotation(Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  replaceable model PipeFriction = 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
                                                                                        extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PartialPipeFriction 
      "Friction model"                                                                                  annotation(choicesAllMatching, Dialog(tab="Advanced", group="Heat transfer and friction loss"));
  replaceable model HeatTransfer = 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      extends 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PartialPipeHeatTransfer
      "Heat transfer model"                                                                       annotation(choicesAllMatching, Dialog(tab="Advanced", group="Heat transfer and friction loss"));
    
  //Display variables
  SI.HeatFlowRate Q_flow_1 "Total heat flow rate of inner pipe";
  SI.HeatFlowRate Q_flow_2 "Total heat flow rate of outer pipe";
    
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe_1(
    redeclare package Medium = Medium_1,
    n=n,
    H_ab=0,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    gravityTerm=false,
    dynamicTerm=false,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    length=length,
    use_wall=true,
    redeclare model Wall = 
        Modelica_Fluid.WorkInProgress.Components.Wall_constProps (
        d_wall=d_wall,
        c_wall=c_wall,
        T_start=T_start_w),
    redeclare HeatTransfer heat,
    redeclare PipeFriction friction,
    initOption=initOption_1,
    use_T_start=use_T_start_1,
    p_start=p_start_1,
    T_start=T_start_1,
    h_start=h_start_1,
    X_start=X_start_1,
    mflow_start=mflow_start_1,
    d_inner=di_1,
    d_outer=da_1)            annotation (extent=[-40,-60; 20,0]);
    
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe_2(
    redeclare package Medium = Medium_2,
    n=n,
    H_ab=0,
    static=static,
    lumped_dp=lumped_dp,
    kineticTerm=kineticTerm,
    gravityTerm=false,
    dynamicTerm=false,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.General,
    length=length,
    redeclare HeatTransfer heat,
    redeclare PipeFriction friction,
    use_T_start=use_T_start_2,
    p_start=p_start_2,
    T_start=T_start_2,
    h_start=h_start_2,
    X_start=X_start_2,
    initOption=initOption_2,
    mflow_start=mflow_start_2,
    P_inner=Modelica.Constants.pi*(da_1 + di_2),
    A_inner=Modelica.Constants.pi/4*(di_2*di_2 - da_1*da_1),
    area_h=Modelica.Constants.pi*da_1*length*K1) 
              annotation (extent=[-40,88; 20,28]);
  annotation (Diagram(Line(points=[-10,36; -10,-8], style(
          color=1,
          rgbcolor={255,0,0},
          thickness=2))), Icon(
      Rectangle(extent=[-100,-26; 100,-30], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,30; 100,26], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95},
          fillPattern=7)),
      Rectangle(extent=[-100,60; 100,30], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,-30; 100,-60], style(
          color=0,
          rgbcolor={0,0,0},
          gradient=2,
          fillColor=70,
          rgbfillColor={0,63,125})),
      Rectangle(extent=[-100,26; 100,-26], style(
          color=69,
          rgbcolor={0,128,255},
          gradient=2,
          fillColor=69,
          rgbfillColor={0,128,255})),
      Text(
        extent=[-100,-60; 100,-100],
        string="%name",
        style(color=3, rgbcolor={0,0,255}))));
  Modelica_Fluid.Interfaces.FluidPort_b port_b1(redeclare package Medium = 
        Medium_1) annotation (extent=[100,-12; 120,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a1(redeclare package Medium = 
        Medium_1) annotation (extent=[-120,-12; -100,8]);
  Modelica_Fluid.Interfaces.FluidPort_a port_a2(redeclare package Medium = 
        Medium_2) annotation (extent=[-120,36; -100,56]);
  Modelica_Fluid.Interfaces.FluidPort_b port_b2(redeclare package Medium = 
        Medium_2) annotation (extent=[100,-56; 120,-36]);
equation 
  Q_flow_1=sum(pipe_1.Qs_flow);
  Q_flow_2=sum(pipe_2.Qs_flow);
  connect(pipe_2.thermalPort, pipe_1.thermalPort);
  connect(pipe_2.port_b, port_b2) annotation (points=[20,58; 60,58; 60,-46; 110,
        -46], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_b, port_b1) annotation (points=[20,-30; 42,-30; 42,-2;
        110,-2], style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_1.port_a, port_a1) annotation (points=[-40.6,-30; -75.3,-30;
          -75.3,-2; -110,-2],
                            style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
  connect(pipe_2.port_a, port_a2) annotation (points=[-40.6,58; -76,58; -76,46;
          -110,46],
                  style(
      color=69,
      rgbcolor={0,127,255},
      thickness=2,
      gradient=2,
      fillColor=42,
      rgbfillColor={213,0,0}));
end HeatExchanger;
  
model Wall_constProps 
    "Pipe wall, assuming ideal 1D-conduction and constant material properties" 
  extends Modelica_Fluid.WorkInProgress.Interfaces.PartialPipeWall;
  parameter SI.Density d_wall "Density of wall material";
  parameter SI.SpecificHeatCapacity c_wall 
      "Specific heat capacity of wall material";
  parameter SI.Temperature T_start "Start value for wall temperature";
  parameter SI.Mass[n] m=ones(n)*(a_outer-a_inner)*length*d_wall/n "Wall mass";
  parameter Modelica_Fluid.Types.Init.Temp initOption;
  SI.Temperature[n] T(start=ones(n)*T_start, stateSelect=StateSelect.prefer) 
      "Wall temperature";
initial equation 
  if initOption==2 then //2: full steady state initialization
    der(T)=zeros(n);
  else
    T=ones(n)*T_start;
  end if;
equation 
  for i in 1:n loop
   c_wall*m[i]*der(T[i]) = thermalPort_a[i].Q_flow + thermalPort_b[i].Q_flow;
  end for;
  //assuming ideal heat conduction perpendicular to fluid flow, conduction in remaining two dimensions is neglected
  thermalPort_a.T=T;
  thermalPort_b.T=T;
end Wall_constProps;
  
  model ShortPipe2 
    "Short pipe with two volumes, wall friction and gravity effect" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid;
    
    extends Modelica_Fluid.Interfaces.PartialTwoPort;
    replaceable package WallFriction = 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.QuadraticTurbulent 
      extends 
      Modelica_Fluid.PressureLosses.Utilities.WallFriction.PartialWallFriction 
      "Characteristic of wall friction"  annotation(choicesAllMatching=true);
    
    parameter SI.Length length "Length of pipe";
    parameter SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    parameter SI.Length height_ab = 0.0 "Height of port_b over port_a" annotation(Evaluate=true);
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
        annotation(Dialog(enable=WallFriction.use_roughness));
    parameter Boolean use_nominal = false 
      "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
        annotation(Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal = 0.01 
      "Nominal dynamic viscosity (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter SI.Density d_nominal = 0.01 
      "Nominal density (for wall friction computation)" annotation(Dialog(enable=use_nominal));
    parameter SI.AbsolutePressure dp_small = 1 
      "Turbulent flow for wall friction if |dp| >= dp_small" 
      annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
    final parameter SI.Volume V = Modelica.Constants.pi*(diameter/2)^2*length;
    
    parameter Types.InitWithGlobalDefault.Temp initVolume1=
              Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option for volume 1" 
      annotation(Dialog(tab = "Initialization"));
    
    parameter Types.InitWithGlobalDefault.Temp initVolume2=
              Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option for volume 2" 
      annotation(Dialog(tab = "Initialization"));
    
    parameter Medium.AbsolutePressure p_start = Medium.p_default 
      "Start value of pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
      "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    annotation (defaultComponentName="pipe",Icon(
        Rectangle(extent=[-100,60; 100,-60],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,34; 100,-36],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Text(
          extent=[-150,-60; 150,-110],
          string="%name",
          style(gradient=2, fillColor=69)),
        Ellipse(extent=[-90,15; -60,-15], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Ellipse(extent=[60,15; 90,-15],   style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))),       Documentation(info="<html>
<p>
Simple pipe model consisting of two volumes, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
</p>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
    PressureLosses.WallFrictionAndGravity frictionAndGravity(
      redeclare package Medium = Medium,
      flowDirection=flowDirection,
      redeclare package WallFriction = WallFriction,
      diameter=diameter,
      roughness=roughness,
      use_nominal=use_nominal,
      eta_nominal=eta_nominal,
      d_nominal=d_nominal,
      from_dp=true,
      dp_small=dp_small,
      show_Re=false,
      length=length,
      height_ab=height_ab) 
                         annotation (extent=[-10,-10; 10,10]);
    Modelica_Fluid.Utilities.PortVolume volume1(
      redeclare package Medium = Medium,
      p_start=p_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initOption=initVolume1) 
      annotation (extent=[-70,-10; -50,10]);
    Modelica_Fluid.Utilities.PortVolume volume2(
      redeclare package Medium = Medium,
      p_start=p_start,
      use_T_start=use_T_start,
      T_start=T_start,
      h_start=h_start,
      X_start=X_start,
      V=V/2,
      initOption=initVolume2) 
      annotation (extent=[50,-10; 70,10]);
  equation 
    connect(volume1.port, port_a) 
      annotation (points=[-60,0; -100,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume1.port, frictionAndGravity.port_a) 
      annotation (points=[-60,0; -10,0], style(color=69, rgbcolor={0,127,255}));
    connect(frictionAndGravity.port_b, volume2.port) 
      annotation (points=[10,0; 60,0], style(color=69, rgbcolor={0,127,255}));
    connect(volume2.port, port_b) 
      annotation (points=[60,0; 100,0], style(color=69, rgbcolor={0,127,255}));
  end ShortPipe2;
  
  model PipeAttachment2 "Equations to attach pipe at tank" 
    import SI = Modelica.SIunits;
      replaceable package Medium = PackageMedium extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    
      Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10,-112; 10,-92],    rotation=90);
     // Real mXi_flow;
      parameter SI.Height pipeHeight "Height of pipe";
      parameter SI.Diameter pipeDiameter "Inner (hydraulic) pipe diameter";
      input SI.AbsolutePressure p_ambient "Tank surface pressure" annotation(Dialog);
      input SI.Height level "Actual tank level" annotation(Dialog);
      input Medium.SpecificEnthalpy h 
      "Actual specific enthalpy of fluid in tank"                 annotation(Dialog);
      input Medium.Density d "Actual specific density of fluid in tank" 
                                                        annotation(Dialog);
      input Medium.MassFraction Xi[Medium.nXi] 
      "Actual mass fractions of fluid in tank"                    annotation(Dialog);
      parameter Real k_small(min=0) = 0 
      "Small regularization range if tank level is below pipe height; k_small = 0 gives ideal switch"
                annotation(Evaluate=true);
      parameter Medium.MassFlowRate m_flow_nominal = 1 
      "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
    
      output Medium.EnthalpyFlowRate H_flow 
      "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate m_flow 
      "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate mXi_flow[Medium.nXi] 
      "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter SI.Area pipeArea = Modelica.Constants.pi*(pipeDiameter/2)^2;
    SI.Length aboveLevel = level - pipeHeight;
    Boolean m_flow_out(start=true) "true= massflow out of tank";
    
  equation 
  m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
    
  if (aboveLevel > 0) then
    port.p = aboveLevel*fluidOptions.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
  else
   if pre(m_flow_out) then
      m_flow = 0;
    else
      port.p = p_ambient;
    end if;
  end if;
  /* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, most examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
 
HTML documentation:
 
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">
 
 
  Real s(start=0) 
    "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
equation 
  m_flow_out = s < 0;
  if aboveLevel > 0 then
     m_flow = m_flow_nominal*s 
      "equation to compute s, which is a dummy in this branch";
     port.p = p_ambient + aboveLevel*fluidOptions.g*d  -
                          smooth(2,if m_flow_out then m_flow^2/(2*d*pipeArea^2) else -k_small*p_ambient*s);
  else
     m_flow = m_flow_nominal*(if m_flow_out then k_small*s else s);
     port.p = p_ambient*(1 + (if m_flow_out then s else k_small*s));
  end if;
*/
    
    /* Handle reverse and zero flow */
    H_flow = port.H_flow;
    m_flow = port.m_flow;
    mXi_flow = port.mXi_flow;
    port.H_flow = semiLinear(port.m_flow, port.h, h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
    annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
</html>"), Icon(Rectangle(extent=[-100,0; 100,-100], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Text(
          extent=[-122,48; 132,6],
          style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1),
          string="%name")));
  end PipeAttachment2;
  
  model OpenTank "Tank with three inlet/outlet-arrays at variable heights" 
    import SI = Modelica.SIunits;
    import Modelica_Fluid.Types;
    replaceable package Medium = PackageMedium 
      extends Modelica.Media.Interfaces.PartialMedium "Medium in the component"
      annotation (choicesAllMatching=true);
    
    parameter SI.Height height "Height of tank";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
    parameter Real k = 4.9 "Heat transfer coefficient from tank to ambient";
    
    parameter Integer n_topPorts = 0 "Number of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; in and out flow of tank)"));
    parameter SI.Height top_heights[n_topPorts]=fill(height,n_topPorts) 
      "Heights of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; in and out flow of tank)"));
    parameter SI.Diameter top_diameters[n_topPorts] 
      "Inner (hydraulic) diameters of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; in and out flow of tank)"));
    
    parameter Integer n_bottomPorts = 0 "Number of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Height bottom_heights[n_bottomPorts]=fill(0.0,n_bottomPorts) 
      "Heights of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Diameter bottom_diameters[n_bottomPorts] 
      "Inner (hydraulic) diameters of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    
    parameter Integer n_sidePorts = 0 "Number of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)"));
    parameter SI.Height side_heights[n_sidePorts] "Heights of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)"));
    parameter SI.Area side_diameters[n_sidePorts] 
      "Inner (hydraulic) diameters of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)"));
    
    parameter Medium.AbsolutePressure p_ambient=fluidOptions.default_p_ambient 
      "Tank surface pressure" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient = fluidOptions.default_T_ambient 
      "Tank surface Temperature" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    
    parameter Types.InitWithGlobalDefault.Temp initOption=
              Types.InitWithGlobalDefault.UseGlobalFluidOption 
      "Initialization option" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start "Start value of tank level" 
      annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon topPorts[n_topPorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,100; 30,108]);
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon bottomPorts[n_bottomPorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,-108; 30,-100],   rotation=90);
    Modelica_Fluid.Interfaces.FluidPort_ArrayIcon sidePorts[n_sidePorts](redeclare 
        package Medium = Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[100,30; 108,-30]);
    
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    
    SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
      "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U "Internal energy of tank volume";
    SI.Mass m "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    Real Q_lost "Wärmeverlust";
    
  encapsulated model PipeAttachment "Equations to attach pipe at tank" 
      import Modelica_Fluid;
      import Modelica;
      import SI = Modelica.SIunits;
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10,-112; 10,-92],    rotation=90);
     // Real mXi_flow;
      parameter SI.Height pipeHeight "Height of pipe";
      parameter SI.Diameter pipeDiameter "Inner (hydraulic) pipe diameter";
      input SI.AbsolutePressure p_ambient "Tank surface pressure" annotation(Dialog);
      input SI.Height level "Actual tank level" annotation(Dialog);
      input Medium.SpecificEnthalpy h 
        "Actual specific enthalpy of fluid in tank"               annotation(Dialog);
      input Medium.Density d "Actual specific density of fluid in tank" 
                                                        annotation(Dialog);
      input Medium.MassFraction Xi[Medium.nXi] 
        "Actual mass fractions of fluid in tank"                  annotation(Dialog);
      parameter Real k_small(min=0) = 0 
        "Small regularization range if tank level is below pipe height; k_small = 0 gives ideal switch"
                annotation(Evaluate=true);
      parameter Medium.MassFlowRate m_flow_nominal = 1 
        "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
      
      output Medium.EnthalpyFlowRate H_flow 
        "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate m_flow 
        "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate mXi_flow[Medium.nXi] 
        "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
    protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
        "Global default options";
    parameter SI.Area pipeArea = Modelica.Constants.pi*(pipeDiameter/2)^2;
    SI.Length aboveLevel = level - pipeHeight;
    Boolean m_flow_out(start=true) "true= massflow out of tank";
      
  equation 
  m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
      
  if (aboveLevel > 0) then
    port.p = aboveLevel*fluidOptions.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
  else
   if pre(m_flow_out) then
      m_flow = 0;
    else
      port.p = p_ambient;
    end if;
  end if;
  /* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, most examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
 
HTML documentation:
 
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">
 
 
  Real s(start=0) 
    "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
equation 
  m_flow_out = s < 0;
  if aboveLevel > 0 then
     m_flow = m_flow_nominal*s 
      "equation to compute s, which is a dummy in this branch";
     port.p = p_ambient + aboveLevel*fluidOptions.g*d  -
                          smooth(2,if m_flow_out then m_flow^2/(2*d*pipeArea^2) else -k_small*p_ambient*s);
  else
     m_flow = m_flow_nominal*(if m_flow_out then k_small*s else s);
     port.p = p_ambient*(1 + (if m_flow_out then s else k_small*s));
  end if;
*/
      
    /* Handle reverse and zero flow */
    H_flow = port.H_flow;
    m_flow = port.m_flow;
    mXi_flow = port.mXi_flow;
    port.H_flow = semiLinear(port.m_flow, port.h, h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
    annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
</html>"), Icon(Rectangle(extent=[-100,0; 100,-100], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Text(
          extent=[-122,48; 132,6],
          style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1),
          string="%name")));
  end PipeAttachment;
    
  protected 
    outer Modelica_Fluid.Components.FluidOptions fluidOptions 
      "Global default options";
    parameter Medium.SpecificEnthalpy h_start = Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start);
    parameter Types.Init.Temp initOption2=
        if initOption == Types.InitWithGlobalDefault.UseGlobalFluidOption then 
             fluidOptions.default_initOption else initOption 
        annotation(Evaluate=true, Hide=true);
    parameter Integer precision = 3 "Precision for tank level in animation" annotation(Hide=false);
    
    Medium.EnthalpyFlowRate H_flow_topPorts[n_topPorts];
    Medium.EnthalpyFlowRate H_flow_bottomPorts[n_bottomPorts];
    Medium.EnthalpyFlowRate H_flow_sidePorts[n_sidePorts];
    
    Medium.MassFlowRate m_flow_topPorts[n_topPorts];
    Medium.MassFlowRate m_flow_bottomPorts[n_bottomPorts];
    Medium.MassFlowRate m_flow_sidePorts[n_sidePorts];
    
    Medium.MassFlowRate mXi_flow_topPorts[n_topPorts,Medium.nXi];
    Medium.MassFlowRate mXi_flow_bottomPorts[n_bottomPorts,Medium.nXi];
    Medium.MassFlowRate mXi_flow_sidePorts[n_sidePorts,Medium.nXi];
    
    OpenTank.PipeAttachment pipeAttachmentTop[n_topPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      H_flow=H_flow_topPorts,
      m_flow=m_flow_topPorts,
      mXi_flow = mXi_flow_topPorts,
      pipeHeight=top_heights,
      pipeDiameter=top_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,80; 20,40]);
    OpenTank.PipeAttachment pipeAttachmentBottom[n_bottomPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      H_flow=H_flow_bottomPorts,
      m_flow=m_flow_bottomPorts,
      mXi_flow = mXi_flow_bottomPorts,
      pipeHeight=bottom_heights,
      pipeDiameter=bottom_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,-80; 20,-40]);
    OpenTank.PipeAttachment pipeAttachmentSide[n_sidePorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      H_flow=H_flow_sidePorts,
      m_flow=m_flow_sidePorts,
      mXi_flow = mXi_flow_sidePorts,
      pipeHeight=side_heights,
      pipeDiameter=side_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[40,-20; 80,20], rotation=90);
    
  equation 
    for i in 1:n_bottomPorts loop
      connect(pipeAttachmentBottom[i].port, bottomPorts[i]) 
        annotation (points=[0,-80.4; 0,-104],
                                            style(color=3, rgbcolor={0,0,255}));
    end for;
    
    for i in 1:n_topPorts loop
       connect(pipeAttachmentTop[i].port, topPorts[i]) 
         annotation (points=[0,80.4; 0,104],
                                           style(color=3, rgbcolor={0,0,255}));
    end for;
    
    for i in 1:n_sidePorts loop
       connect(pipeAttachmentSide[i].port, sidePorts[i]) 
         annotation (points=[80.4,-1.2491e-015; 104,0],
                                             style(color=3, rgbcolor={0,0,255}));
    end for;
    
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    Q_lost = - k*2*sqrt(Modelica.Constants.pi*area)*level*(medium.T - T_ambient);
    
    // Mass balances
    der(m) = sum(m_flow_bottomPorts) + sum(m_flow_sidePorts) + sum(m_flow_topPorts);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(mXi_flow_bottomPorts[i,:]) +
                    sum(mXi_flow_sidePorts[i,:]) +
                    sum(mXi_flow_topPorts[i,:]);
    end for;
    
    // Energy balance
    if Medium.singleState then
      der(U) = sum(H_flow_bottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) + Q_lost 
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      der(U) = sum(H_flow_bottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) - p_ambient*der(V) + Q_lost;
    end if;
    
     assert(level < height,   " 
    Tank ist überfüllt.
    ");
    
  initial equation 
    if initOption2 == Types.Init.NoInit then
      // no initial equations
    elseif initOption2 == Types.Init.InitialValues then
      level = level_start;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initOption2 == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initOption2 == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
    annotation (defaultComponentName="tank",
      Icon(
        Rectangle(extent=[-100,100; 100,0],     style(color=7, fillColor=7)),
        Rectangle(extent=DynamicSelect([-100,-100; 100,10],
                                       [-100, -100; 100, (-100 + 200*level/height)]),
            style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=71,
            rgbfillColor={85,170,255},
            fillPattern=1)),
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,90; 95,60],
          string="%name",
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,50; 95,35],
          style(color=0),
          string="%level_start",
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3))),
        Text(
          extent=[-95,30; 95,5],
          style(color=0),
          string=DynamicSelect(" ",realString(level,1,integer(precision)))),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3)),
        Text(
          extent=[-40,124; -22,112],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[111,37; 129,25],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[-38,-112; -20,-124],
          string="1",
          style(color=10, rgbcolor={135,135,135}))),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. The top part is open to the environment at the fixed pressure <tt>p_ambient</tt>. Heat transfer to the environment and to the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, assumed to have uniform temperature and mass fractions.
<p>The geometry of the tank is specified by the following parameters: <tt>V0</tt> is the volume of the liquid when the level is at the zero reference; <tt>area</tt> is the cross-sectional area of the tank; <tt>H0</tt> is the height of the zero-reference level plane over the port connector. It is thus possible to model rounded-bottom tanks, as long as they have a cylindrical shape in the range of operating levels.
<p>The tank can be initialized with the following options:
<ul>
<li>NoInit: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
 
 
 
</HTML>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
  end OpenTank;
  
end Components;
