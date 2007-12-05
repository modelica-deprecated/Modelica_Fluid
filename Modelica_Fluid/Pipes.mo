within Modelica_Fluid;
package Pipes "Lumped, distributed and thermal pipe components" 
    extends Modelica_Fluid.Icons.VariantLibrary;
  
  model DistributedPipe_obsolete "Distributed pipe model" 
    
    extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow_obsolete(
      Qs_flow=heatTransfer.Q_flow,
      Ws_flow=zeros(n),
      ms_flow=zeros(n),
      msXi_flow=zeros(n, Medium.nXi),
      Vi=ones(n)*V/n);
    parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));
    final parameter SI.Volume V=area*length;
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
    
  //Pressure Drop
    parameter Boolean kineticTerm=false 
      " = true, include kinetic term in momentum balance" annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                  extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
    parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
    parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                   annotation(Dialog(tab="General", group="Pressure loss"));
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
    parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
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
    inner Medium.ThermodynamicState[n] state;
    replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) extends 
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) "Convective heat transfer" 
              annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;
      20,20]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
  annotation (extent=[-10,44; 10,64]);
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  protected 
    SI.Pressure[n+1] dp_stat;
  //  Real DI_flow[n+1] "Momentum flow (if kineticTerm=true)";
    SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
        use_eta_nominal then eta_nominal else eta[1]);
    SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else eta[n]);
    SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
        Medium.dynamicViscosity(medium));
    
  annotation (
    Icon(
      Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
    Diagram,
    Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. 
</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    
  //Pressure difference in momentum balance
  dp = dp_stat;
    // dp = DI_flow + dp_stat;
    // state=medium.state;
    // if kineticTerm then
    //   for i in 2:n loop
    //     DI_flow[i] = Modelica_Fluid.Utilities.regSquare(m_flow[i], delta=max(
    //       0.0001, 0.01*mflow_start))/area*(1/d[i - 1] - 1/d[i]);
    //   end for;
    //   DI_flow[1] = Modelica_Fluid.Utilities.regSquare(m_flow[1], delta=max(0.0001,
    //     0.01*mflow_start))/area*(1/d_a - 1/d[1]);
    //   DI_flow[n+1] = Modelica_Fluid.Utilities.regSquare(m_flow[n+1], delta=max(
    //     0.0001, 0.01*mflow_start))/area*(1/d[n] - 1/d_b);
    // else
    //   DI_flow = fill(0, n+1);
    // end if;
    
    //Pressure drop and gravity
  if from_dp and not WallFriction.dp_is_zero then
     m_flow[1] = WallFriction.massFlowRate_dp(
        dp_stat[1] - height_ab/2/n*ambient.g*(d_a + d[1])/2,
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n/2,
        d_h,
        roughness,
        dp_small);
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp(
          dp_stat[i] - height_ab/n*ambient.g*(d[i - 1] + d[i])/2,
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          dp_small);
      end for;
      m_flow[n + 1] = WallFriction.massFlowRate_dp(
        dp_stat[n+1] - height_ab/n/2*ambient.g*(d[n] + d_b)/2,
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        dp_small);
      
  else
      
      dp_stat[1] = WallFriction.pressureLoss_m_flow(
        m_flow[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/2/n*ambient.g*(d_a + d[1])/2;
      for i in 2:n loop
        dp_stat[i] = WallFriction.pressureLoss_m_flow(
          m_flow[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
      end for;
      dp_stat[n+1] = WallFriction.pressureLoss_m_flow(
        m_flow[n+1],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/n/2*ambient.g*(d[n] + d_b)/2;
      
  end if;
  connect(thermalPort, heatTransfer.thermalPort) 
    annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_obsolete;
  
 model DistributedPipeLumpedPressure "Distributed pipe model" 
    
   extends 
      Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlowLumpedPressure(
     Qs_flow=heatTransfer.Q_flow,
     Ws_flow=zeros(n),
     ms_flow=zeros(n),
     msXi_flow=zeros(n, Medium.nXi),
     Vi=ones(n)*V/n);
    
    //Volume size
   final parameter SI.Volume V=area*length;
 //Heat transfer area
   parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                              annotation(Dialog(tab="General", group="Heat transfer"));
    //Pressure Drop
   replaceable package WallFriction = 
       Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
     extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction" annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
   parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                  annotation(Dialog(tab="General", group="Pressure loss"));
   parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
   parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
   parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
    
   parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
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
       diameter) if                                                                                              show_Re 
      "Reynolds number of pipe flow";
   inner Medium.ThermodynamicState[n] state;
   replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
      heatTransfer(
     redeclare final package Medium = Medium,
     final n=n,
     final d_h=d_h,
     final A_h=area_h,
     final A_cross=area,
     final L=length,
     T=medium.T) extends 
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
     redeclare final package Medium = Medium,
     final n=n,
     final d_h=d_h,
     final A_h=area_h,
     final A_cross=area,
     final L=length,
     T=medium.T) "Convective heat transfer" 
             annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;
     20,20]);
   Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
 annotation (extent=[-10,44; 10,64]);
   outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  protected 
   SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
       use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[1].state) else (if m_flow[1]>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p, port_a.h, port_a.Xi)) else Medium.dynamicViscosity(medium[1].state))));
   SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then Medium.dynamicViscosity(medium[n].state) else (if m_flow[1]<0 then Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h, port_b.Xi)) else Medium.dynamicViscosity(medium[n].state))));
    
   annotation (
     Icon(
       Rectangle(extent=[-100,44; 100,40], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=10,
           rgbfillColor={95,95,95})),
       Rectangle(extent=[-100,-40; 100,-44], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=10,
           rgbfillColor={95,95,95})),
       Ellipse(extent=[-72,10; -52,-10], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=0,
           rgbfillColor={0,0,0})),
       Ellipse(extent=[-30,10; -10,-10], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=0,
           rgbfillColor={0,0,0})),
       Ellipse(extent=[10,10; 30,-10], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=0,
           rgbfillColor={0,0,0})),
       Ellipse(extent=[50,10; 70,-10], style(
           color=0,
           rgbcolor={0,0,0},
           fillColor=0,
           rgbfillColor={0,0,0})),
       Text(
         extent=[-143,-42; 157,-92],
         string="%name",
         style(gradient=2, fillColor=69))),
     Diagram,
     Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow_pLumped\">PartialDistributedFlow_pLumped</a>. Source terms in mass and energy balances are set to zero. The total volume is a paramter. The number of momentum balances is reduced to two, one on each side of the hydraulic state, which corresponds to a constant pressure along the entire pipe with pressure drop and gravitational forces lumped at the ports.<The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. .
</html>",  revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    
 equation 
   state=medium.state;
    
 if from_dp and not WallFriction.dp_is_zero then
     m_flow[1] = WallFriction.massFlowRate_dp(
       dp[1] - ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*ambient.g*(d_a + d_b)
         /2,
       d_a,
       d_b,
       eta_a,
       eta_b,
       ((integer(n/2) + 1)*2 - 1)/(2*n)*length,
       d_h,
       roughness,
       dp_small);
     m_flow[n + 1] = WallFriction.massFlowRate_dp(
       dp[2] - (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*ambient.g*(d_a
          + d_b)/2,
       d_a,
       d_b,
       eta_a,
       eta_b,
       (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
       d_h,
       roughness,
       dp_small);
 else
   dp[1] = (WallFriction.pressureLoss_m_flow(
       m_flow[1],
       d_a,
       d_b,
       eta_a,
       eta_b,
       ((integer(n/2) + 1)*2 - 1)/(2*n)*length,
       d_h,
       roughness,
       m_flow_small) + ((integer(n/2) + 1)*2 - 1)/(2*n)*height_ab*ambient.g*(
       d_a + d_b)/2);
     dp[2] = (WallFriction.pressureLoss_m_flow(
       m_flow[1],
       d_a,
       d_b,
       eta_a,
       eta_b,
       (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*length,
       d_h,
       roughness,
       m_flow_small) + (2*n - (integer(n/2) + 1)*2 + 1)/(2*n)*height_ab*
       ambient.g*(d_a + d_b)/2);
 end if;
    
   connect(thermalPort, heatTransfer.thermalPort) 
     annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
 end DistributedPipeLumpedPressure;
  
  model DistributedPipe_a_v_b "Distributed pipe model" 
    extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
      Qs_flow=heatTransfer.Q_flow,
      Ws_flow=zeros(n),
      ms_flow=zeros(n),
      msXi_flow=zeros(n, Medium.nXi),
      Vi=ones(n)*V/n,
      final modelStructure=Modelica_Fluid.Types.ModelStructure.a_v_b,
      redeclare final Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare 
          package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
      redeclare final Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare 
          package Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));
    
    parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));
    final parameter SI.Volume V=area*length;
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
    
  //Pressure Drop
    replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                  extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
    parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
    parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                   annotation(Dialog(tab="General", group="Pressure loss"));
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
    parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
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
    inner Medium.ThermodynamicState[n] state;
    replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) extends 
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) "Convective heat transfer" 
              annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;
      20,20]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
  annotation (extent=[-10,44; 10,64]);
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  protected 
    SI.Pressure[n+1] dp_stat;
    SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
        use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else (if m_flow[1]>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p, port_a.h, port_a.Xi)) else eta[1])));
    SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else (if m_flow[1]<0 then Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h, port_b.Xi)) else eta[n])));
    SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
        Medium.dynamicViscosity(medium.state));
    
  annotation (
    Icon(
      Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
    Diagram,
    Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>

</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    
  //Pressure difference in momentum balance
  dp=dp_stat;
    state=medium.state;
    
    //Pressure drop and gravity
  if from_dp and not WallFriction.dp_is_zero then
      m_flow[1] = WallFriction.massFlowRate_dp(
        dp_stat[1] - height_ab/2/n*ambient.g*d[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n/2,
        d_h,
        roughness,
        dp_small);
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp(
          dp_stat[i] - height_ab/n*ambient.g*(d[i-1] + d[i])/2,
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          dp_small);
      end for;
      m_flow[n + 1] = WallFriction.massFlowRate_dp(
        dp_stat[n+1] - height_ab/n/2*ambient.g*d[n],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        dp_small);
      
  else
       dp_stat[1] = WallFriction.pressureLoss_m_flow(
          m_flow[1],
          d_a,
          d[1],
          eta_a,
          eta[1],
          length/n/2,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*d[1]/2;
      for i in 2:n loop
        dp_stat[i] = WallFriction.pressureLoss_m_flow(
          m_flow[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
      end for;
      dp_stat[n] = WallFriction.pressureLoss_m_flow(
        m_flow[n+1],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n/2,
        d_h,
        roughness,
        m_flow_small) + height_ab/n/2*ambient.g*d[n]/2;
      
  end if;
  connect(thermalPort, heatTransfer.thermalPort) 
    annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_a_v_b;
  
  model DistributedPipe_a_vb "Distributed pipe model" 
    
    extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
      Qs_flow=heatTransfer.Q_flow,
      Ws_flow=zeros(n),
      ms_flow=zeros(n),
      msXi_flow=zeros(n, Medium.nXi),
      Vi=ones(n)*V/n,
      final modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
      redeclare final Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare 
          package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
      redeclare final Modelica_Fluid.Interfaces.FluidStatePort_b port_b(redeclare 
          package Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));
    parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));
    final parameter SI.Volume V=area*length;
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
    
  //Pressure Drop
   replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                  extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
    parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
    parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                   annotation(Dialog(tab="General", group="Pressure loss"));
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
    parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
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
    inner Medium.ThermodynamicState[n] state;
    replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) extends 
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) "Convective heat transfer" 
              annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;
      20,20]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
  annotation (extent=[-10,44; 10,64]);
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  protected 
    SI.Pressure[n] dp_stat;
     SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
        use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else (if m_flow[1]>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p, port_a.h, port_a.Xi)) else eta[1])));
    SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else (if m_flow[1]<0 then Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h, port_b.Xi)) else eta[n])));
   SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
        Medium.dynamicViscosity(medium.state));
    
  annotation (
    Icon(
      Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
    Diagram,
    Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. 
</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    
  //Pressure difference in momentum balance
  dp=dp_stat;
    state=medium.state;
    
    //Pressure drop and gravity
  if from_dp and not WallFriction.dp_is_zero then
      m_flow[1] = WallFriction.massFlowRate_dp(
        dp_stat[1] - height_ab/n*ambient.g*d[1],
        d_a,
        d[1],
        eta_a,
        eta[1],
        length/n,
        d_h,
        roughness,
        dp_small);
      for i in 2:n loop
        m_flow[i] = WallFriction.massFlowRate_dp(
          dp_stat[i] - height_ab/n*ambient.g*(d[i-1] + d[i])/2,
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          dp_small);
      end for;
  else
       dp_stat[1] = WallFriction.pressureLoss_m_flow(
          m_flow[1],
          d_a,
          d[1],
          eta_a,
          eta[1],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*d[1];
      for i in 2:n loop
        dp_stat[i] = WallFriction.pressureLoss_m_flow(
          m_flow[i],
          d[i - 1],
          d[i],
          eta[i - 1],
          eta[i],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
      end for;
  end if;
  connect(thermalPort, heatTransfer.thermalPort) 
    annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_a_vb;
  
  model DistributedPipe_av_b "Distributed pipe model" 
    
    extends Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow(
      Qs_flow=heatTransfer.Q_flow,
      Ws_flow=zeros(n),
      ms_flow=zeros(n),
      msXi_flow=zeros(n, Medium.nXi),
      Vi=ones(n)*V/n,
      final modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
      redeclare final Modelica_Fluid.Interfaces.FluidStatePort_a port_a(redeclare 
          package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -Modelica.Constants.inf else  0)),
      redeclare final Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare 
          package Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +Modelica.Constants.inf else 0)));
    parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                               annotation(Dialog(tab="General", group="Heat transfer"));
    final parameter SI.Volume V=area*length;
    parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                     annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
    
  //Pressure Drop
   replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
                                                                                  extends 
      Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
      "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
    parameter Boolean use_d_nominal=false 
      "= true, if d_nominal is used, otherwise computed from medium"                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter SI.Density d_nominal=0.01 
      "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
    parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                                   annotation(Dialog(tab="General", group="Pressure loss"));
    parameter SI.Length roughness(min=0) = 2.5e-5 
      "Absolute roughness of pipe (default = smooth steel pipe)" 
      annotation(Dialog(tab="General", group="Pressure loss",enable=WallFriction.use_roughness));
    parameter Boolean use_eta_nominal=false 
      "= true, if eta_nominal is used, otherwise computed from medium"                            annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
    parameter SI.DynamicViscosity eta_nominal=0.01 
      "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
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
    inner Medium.ThermodynamicState[n] state;
    replaceable Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha
      heatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) extends 
      Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
      redeclare final package Medium = Medium,
      final n=n,
      final d_h=d_h,
      final A_h=area_h,
      final A_cross=area,
      final L=length,
      T=medium.T) "Convective heat transfer" 
              annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;
      20,20]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
      "Thermal port" 
  annotation (extent=[-10,44; 10,64]);
    outer Modelica_Fluid.Ambient ambient "Ambient conditions";
    
  protected 
    SI.Pressure[n] dp_stat;
     SI.DynamicViscosity eta_a=if not WallFriction.use_eta then 1.e-10 else (if 
        use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[1] else (if m_flow[1]>=0 then Medium.dynamicViscosity(Medium.setState_phX(port_a.p, port_a.h, port_a.Xi)) else eta[1])));
    SI.DynamicViscosity eta_b=if not WallFriction.use_eta then 1.e-10 else (if use_eta_nominal then eta_nominal else (if use_approxPortProperties then eta[n] else (if m_flow[1]<0 then Medium.dynamicViscosity(Medium.setState_phX(port_b.p, port_b.h, port_b.Xi)) else eta[n])));
    SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else (if use_eta_nominal then fill(eta_nominal, n) else 
        Medium.dynamicViscosity(medium.state));
    
  annotation (
    Icon(
      Rectangle(extent=[-100,44; 100,40], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Rectangle(extent=[-100,-40; 100,-44], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=10,
          rgbfillColor={95,95,95})),
      Ellipse(extent=[-72,10; -52,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[-30,10; -10,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[10,10; 30,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Ellipse(extent=[50,10; 70,-10], style(
          color=0,
          rgbcolor={0,0,0},
          fillColor=0,
          rgbfillColor={0,0,0})),
      Text(
        extent=[-143,-42; 157,-92],
        string="%name",
        style(gradient=2, fillColor=69))),
    Diagram,
    Documentation(info="<html>
Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>. 
</html>",
        revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  equation 
    
  //Pressure difference in momentum balance
  dp=dp_stat;
    state=medium.state;
    
    //Pressure drop and gravity
  if from_dp and not WallFriction.dp_is_zero then
      for i in 1:n-1 loop
        m_flow[i+1] = WallFriction.massFlowRate_dp(
          dp_stat[i] - height_ab/n*ambient.g*(d[i] + d[i+1])/2,
          d[i],
          d[i+1],
          eta[i],
          eta[i+1],
          length/n,
          d_h,
          roughness,
          dp_small);
      end for;
      m_flow[n+1] = WallFriction.massFlowRate_dp(
        dp_stat[n] - height_ab/n*ambient.g*d[n],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n,
        d_h,
        roughness,
        dp_small);
      
  else
      for i in 1:n-1 loop
        dp_stat[i] = WallFriction.pressureLoss_m_flow(
          m_flow[i+1],
          d[i],
          d[i+1],
          eta[i],
          eta[i+1],
          length/n,
          d_h,
          roughness,
          m_flow_small) + height_ab/n*ambient.g*(d[i] + d[i+1])/2;
      end for;
      dp_stat[n] = WallFriction.pressureLoss_m_flow(
        m_flow[n+1],
        d[n],
        d_b,
        eta[n],
        eta_b,
        length/n,
        d_h,
        roughness,
        m_flow_small) + height_ab/n/2*ambient.g*d[n];
  end if;
  connect(thermalPort, heatTransfer.thermalPort) 
    annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
  end DistributedPipe_av_b;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    
  partial model PartialDistributedFlow_obsolete 
      import Modelica_Fluid.Types;
      import Modelica.Constants.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
        annotation (choicesAllMatching=true);
      
  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
      
  //Advanced model options
    parameter Boolean allowFlowReversal=true 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                       annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
    parameter Boolean static=false 
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
      
  //Initialization
      parameter Types.Init.Temp initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start" 
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
      parameter Medium.MassFlowRate mflow_start 
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      
  //Geometry
    parameter Boolean isCircular=true 
        "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0 
        "Inner perimeter"                                                                                     annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0 
        "Inner cross section area"          annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";
      
  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Fluid mass";
    SI.Mass[n,Medium.nXi] mXi "Substance mass";
      
  //Flow quantities
    inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -inf else 
                0, each start=mflow_start, each fixed=false) 
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "Velocity at volume boundaries (not used in balances)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow 
        "Enthalpy flow rates of fluid across segment boundaries";
   parameter Boolean use_d_nominal=false 
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal=0.01 
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
   SI.Pressure[n+1] dp(start=dp0) "pressure difference across staggered grid";
      
  //Fluid ports
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
        "Fluid inlet port" 
                         annotation (extent=[-110,-10; -90,10]);
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
      
     annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69))),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). Half a momentum balance is used on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows taking into account friction losses with respect to the actual geometry of the connection point.</p>
 

 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      
    //Source terms, have to be set in inheriting class (to zero if not used)
    protected 
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
        "Independent mass flow rates, source or sink";
    SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
      
    final parameter SI.Pressure[n+1] dp0={dp_start/(if i
         > 1 and i < n+1 then n else 2*n) for i in 1:n+1} 
        "pressure difference start values";
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
   //SI.Density d_a=if use_d_nominal then d_nominal else d[1];//approximation
    //SI.Density d_b=if use_d_nominal then d_nominal else d[n];//approximation
    SI.Density d_a=if m_flow[1]>=0 then Medium.density_phX(port_a.p, port_a.h, port_a.Xi) else d[1];
    SI.Density d_b=if m_flow[n+1]>=0 then d[n] else Medium.density_phX(port_b.p, port_b.h, port_b.Xi);
      
  equation 
    // Boundary conditions
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
    port_a.mC_flow = semiLinear(port_a.m_flow, port_a.C, port_b.C);
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
      
    // Distributed flow quantities, upwind discretization
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
      
    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;
      
    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
     else
    //dynamic mass and energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    end if;
    zeros(Medium.nC)=port_a.mC_flow + port_b.mC_flow;
    for i in 1:n loop
      assert((allowFlowReversal and not static) or (m_flow[i] >= 0), "Flow reversal not allowed in distributed pipe");
    end for;
      
  //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
  //must be specified in extending class
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      dp[n+1]=medium[n].p-port_b.p;
      
  initial equation 
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        if not (Medium.singleState) then
          der(medium.p) = zeros(n);
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
           medium.p=p_start;
        end if;
      elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
          der(medium.p) = zeros(n);
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlow_obsolete;
    
  partial model PartialDistributedFlowLumpedPressure 
      "Base class for 1D fluid flow with the number of momentum balances reduced to 2" 
    import Modelica_Fluid.Types;
    import Modelica.Constants.*;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
        annotation (choicesAllMatching=true);
      
  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
    final parameter Integer np=2 "Number of momentum balances"                            annotation(Dialog(tab="Advanced"),Evaluate=true);
    final parameter Integer nl=integer(n/2)+1 
        "Number of control volume that contains single state"                 annotation(Evaluate=true);
      
  //Advanced model options
    parameter Boolean allowFlowReversal=true 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                       annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
    parameter Boolean static=false 
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
      
    parameter Boolean use_approxPortProperties=false 
        "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                       annotation(Dialog(tab="Advanced", group="Momentum balances"),Evaluate=true);
      
  //Initialization
      parameter Types.Init.Temp initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
        "Port a pressure start value" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
        "Port b pressure start value" 
      annotation(Dialog(tab = "Initialization"));
      final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Pressure start values";
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
      parameter Medium.MassFlowRate mflow_start 
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      
  //Geometry parameters
    parameter Boolean isCircular=true 
        "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0 
        "Inner perimeter"                                                                                     annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0 
        "Inner cross section area"          annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    SI.Volume[n] Vi "Discretized volume, to be determined in inheriting class ";
      
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
   parameter Boolean use_d_nominal=false 
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal=0.01 
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
      
    SI.Pressure[np] dp(start=dp0) "Pressure difference across staggered grid";
      
  //Fluid ports
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
        "Fluid inlet port" 
                         annotation (extent=[-110,-10; -90,10]);
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
      
     annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69))),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). For this model only two momentum balances are formed on each side of a single pressure state (roughly half way along the flowpath). This assumes a constant pressure level for all medium models in the pipe. The total pressure drop (or rise) is split to be located on each end of the component. Connecting two pipes results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) per momentum balance is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
 
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      
    //Source terms, have to be set in inheriting class
    protected 
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
        "Independent mass flow rates, source or sink";
    SI.HeatFlowRate[n] Qs_flow 
        "Heat and/or enthalpy flow rate, perpendicular to flow direction";
    SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
      
    final parameter SI.Pressure[np] dp0=fill(dp_start/np,np);
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    SI.Density d_a=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[1] else (if m_flow[1]>=0 then Medium.density_phX(port_a.p, port_a.h, port_a.Xi) else d[1]));
    SI.Density d_b=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[n] else (if m_flow[n+1]>=0 then d[n] else Medium.density_phX(port_b.p, port_b.h, port_b.Xi)));
      
  equation 
    // Boundary conditions
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
    port_a.mC_flow = semiLinear(port_a.m_flow, port_a.C, port_b.C);
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
      
    // Distributed flow quantities
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
      
    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;
      
    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i] + Ws_flow[i];
      end for;
    else
     //dynamic mass and energy balances, n "thermal" states, 1 pressure state 
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i] + Ws_flow[i];
      end for;
      end if;
      zeros(Medium.nC) = port_a.mC_flow + port_b.mC_flow;
      for i in 1:n loop
        assert((allowFlowReversal and not static) or (m_flow[i] >= 0),
          "Flow reversal not allowed in distributed pipe");
      end for;
      
  //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
  //two momentum balances, one on each side of pressure state, dp must be supplied in extending class
      dp[1] = port_a.p - medium[nl].p;
      dp[2] = medium[nl].p - port_b.p;
      if n == 2 then
        medium[2].p = medium[1].p;
      elseif n > 2 then
        medium[1:nl - 1].p = ones(nl - 1)*medium[nl].p;
        medium[nl + 1:n].p = ones(n - nl)*medium[nl].p;
      end if;
      
  initial equation 
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        if not Medium.singleState then
          der(medium[nl].p) = 0;
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
        medium[nl].p=p_start[nl];
        end if;
     elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
        der(medium[nl].p) = 0;
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlowLumpedPressure;
    
  partial model PartialDistributedFlow 
    import Modelica_Fluid.Types;
      
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Fluid medium model" 
        annotation (choicesAllMatching=true);
      
  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
      
  //Advanced model options
    parameter Boolean allowFlowReversal=true 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                       annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
    parameter Boolean static=false 
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
      
    parameter Types.ModelStructure.Temp modelStructure=Types.ModelStructure.a_v_b 
        "Determines whether flow model between volume and port is present"                                                                           annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
      
    parameter Boolean use_approxPortProperties=false 
        "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                       annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
      
  //Initialization
      parameter Types.Init.Temp initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
        "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
      final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start,
          p_b_start + (p_a_start - p_b_start)/n,
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
      parameter Medium.MassFlowRate mflow_start 
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      
  //Geometry
    parameter Boolean isCircular=true 
        "= true if cross sectional area is circular" 
      annotation (Evaluate, Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=isCircular));
    parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0 
        "Inner perimeter"                                                                                     annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0 
        "Inner cross section area"          annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
    SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";
      
  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Fluid mass";
    SI.Mass[n,Medium.nXi] mXi "Substance mass";
      
  //Flow quantities
    inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -Modelica.Constants.inf else 
                0, each start=mflow_start, each fixed=false) 
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "Velocity at volume boundaries (not used in balances)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow 
        "Enthalpy flow rates of fluid across segment boundaries";
   parameter Boolean use_d_nominal=false 
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal=0.01 
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
   SI.Pressure[np] dp(start=dp0) "pressure difference across staggered grid";
      
  //Fluid ports
    replaceable Modelica_Fluid.Interfaces.FluidPort port_a "Fluid inlet port" 
                         annotation (extent=[-110,-10; -90,10]);
    replaceable Modelica_Fluid.Interfaces.FluidPort port_b "Fluid outlet port" 
                          annotation (extent=[90,-10; 110,10]);
    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));
      
     annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69))),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). Half a momentum balance is used on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
 
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      
    //Source terms, have to be set in inheriting class (to zero if not used)
    protected 
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
        "Independent mass flow rates, source or sink";
    SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
    final parameter Integer np=if modelStructure==Types.ModelStructure.avb then n-1 else if (modelStructure==Types.ModelStructure.a_vb or modelStructure==Types.ModelStructure.av_b) then n else n+1;
    final parameter SI.Pressure[np] dp0=fill(dp_start,np) 
        "pressure difference start values";
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    SI.Density d_a=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[1] else (if m_flow[1]>=0 then Medium.density_phX(port_a.p, port_a.h, port_a.Xi) else d[1]));
    SI.Density d_b=if use_d_nominal then d_nominal else (if use_approxPortProperties then d[n] else (if m_flow[n+1]>=0 then d[n] else Medium.density_phX(port_b.p, port_b.h, port_b.Xi)));
      
  equation 
    // Boundary conditions
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
    port_a.mC_flow = semiLinear(port_a.m_flow, port_a.C, port_b.C);
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
      
    // Distributed flow quantities, upwind discretization
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
      
    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;
      
    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
     else
    //dynamic mass and energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    end if;
    zeros(Medium.nC)=port_a.mC_flow + port_b.mC_flow;
    for i in 1:n loop
      assert((allowFlowReversal and not static) or (m_flow[i] >= 0), "Flow reversal not allowed in distributed pipe");
    end for;
      
  //Simplified Momentum Balance, dp contains contributions from gravitational and friction effects
  //must be specified in extending class
  if modelStructure==Types.ModelStructure.av_b then
      port_a.p = medium[1].p;
      for i in 1:n-1 loop
        dp[i]=medium[i].p-medium[i+1].p;
      end for;
      dp[n]=medium[n].p-port_b.p;
  elseif modelStructure==Types.ModelStructure.a_vb then
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      port_b.p=medium[n].p;
  elseif modelStructure==Types.ModelStructure.avb then
      port_a.p=medium[1].p;
      for i in 1:n-1 loop
        dp[i]=medium[i].p-medium[i+1].p;
      end for;
      port_b.p=medium[n].p;
  else //ModelStructure.a_v_b
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      dp[n+1]=medium[n].p-port_b.p;
  end if;
      
  initial equation 
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        if not (Medium.singleState) then
          der(medium.p) = zeros(n);
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
           medium.p=p_start;
        end if;
      elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not Medium.singleState then
          der(medium.p) = zeros(n);
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlow;
    
    package CharacteristicNumbers 
      function ReynoldsNumber 
        input SI.MassFlowRate m_flow "Mass flow rate";
        input SI.Length d_ch "Characteristic length (hyd. diam. in pipes)";
        input SI.Area A "Cross sectional area";
        input SI.DynamicViscosity eta "Dynamic viscosity";
        output SI.ReynoldsNumber Re "Reynolds number";
        annotation (Documentation(info="Calculate Re-Number; Re = mdot*Dhyd/A/eta"),
             Icon);
      algorithm 
        Re := abs(m_flow)*d_ch/A/eta;
      end ReynoldsNumber;
      
      function NusseltNumber 
        input SI.CoefficientOfHeatTransfer alpha "Coefficient of heat transfer";
        input SI.Length d_ch "Characteristic length";
        input SI.ThermalConductivity lambda "Thermal conductivity";
        output SI.NusseltNumber Nu "Nusselt number";
        annotation (Documentation(info="Nusselt number Nu = alpha*d_ch/lambda"));
      algorithm 
        Nu := alpha*d_ch/lambda;
      end NusseltNumber;
    end CharacteristicNumbers;
    
  package HeatTransfer 
    partial model PartialPipeHeatTransfer 
        "base class for any pipe heat transfer correlation" 
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
      SI.HeatFlowRate[n] Q_flow "Heat flow rates";
      parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area A_cross "Cross flow area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length L "Total pipe length" annotation(Dialog(tab="No input", enable=false));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
          "Thermal port" 
        annotation (extent=[-20,60; 20,80]);
      input SI.Temperature[n] T;
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
Base class for heat transfer models that can be used in distributed pipe models.
</html>"));
    end PartialPipeHeatTransfer;
      
    partial model PipeHT_Nu 
        "base class for pipe heat transfer correlation in terms of Nusselt numberheat transfer in a circular pipe for laminar and turbulent one-phase flow" 
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=100;
      outer input Medium.ThermodynamicState[n] state;
      outer input SI.MassFlowRate[n+1] m_flow;
      SI.CoefficientOfHeatTransfer[n] alpha(each start=alpha0) 
          "CoefficientOfHeatTransfer";
      Real[n] Re "Reynolds number";
      Real[n] Pr "Prandtl number";
      Real[n] Nu "Nusselt number";
      SI.DynamicViscosity[n] eta "Dynamic viscosity";
      SI.ThermalConductivity[n] lambda "Thermal conductivity";
    equation 
      eta=Medium.dynamicViscosity(state);
      lambda=Medium.thermalConductivity(state);
      Pr = Medium.prandtlNumber(state);
      Re = CharacteristicNumbers.ReynoldsNumber(m_flow[1:n], fill(d_h,n), fill(A_cross, n), eta);
      Nu = CharacteristicNumbers.NusseltNumber(alpha, fill(d_h,n), lambda);
      thermalPort.Q_flow=Q_flow;
      for i in 1:n loop
        thermalPort[i].Q_flow=alpha[i]*A_h/n*(thermalPort[i].T-T[i]);
      end for;
        annotation (Documentation(info="<html>
Base class for heat transfer models that are expressed in terms of the Nusselt number and which can be used in distributed pipe models.
</html>"));
    end PipeHT_Nu;
      
    model PipeHT_constAlpha 
      extends PartialPipeHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0=200;
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a distributed pipe models.
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
      
    model PipeHT_LamTurb_local 
        "laminar and turbulent forced convection in pipes, local coeff." 
      extends PipeHT_Nu;
      protected 
      Real[n] Nu_turb "Nusselt number for turbulent flow";
      Real[n] Nu_lam "Nusselt number for laminar flow";
      Real Nu_1;
      Real[n] Nu_2;
      Real[n] Xi;
      parameter SI.Length dx=L/n;
    equation 
      Nu_1=3.66;
      for i in 1:n loop
       Nu_turb[i]=smooth(0,(Xi[i]/8)*abs(Re[i])*Pr[i]/(1+12.7*(Xi[i]/8)^0.5*(Pr[i]^(2/3)-1))*(1+1/3*(d_h/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(2/3)));
       Xi[i]=(1.8*Modelica.Math.log10(max(1e-10,Re[i]))-1.5)^(-2);
       Nu_lam[i]=(Nu_1^3+0.7^3+(Nu_2[i]-0.7)^3)^(1/3);
       Nu_2[i]=smooth(0,1.077*(abs(Re[i])*Pr[i]*d_h/dx/(if m_flow[i]>=0 then (i-0.5) else (n-i+0.5)))^(1/3));
       Nu[i]=spliceFunction(Nu_turb[i], Nu_lam[i], Re[i]-6150, 3850);
    end for;
      annotation (Documentation(info="<html>
Heat transfer model for laminar and turbulent flow in pipes. Range of validity:
<ul>
<li>fully developed pipe flow</li>
<li>forced convection</li>
<li>one phase Newtonian fluid</li>
<li>(spatial) constant wall temperature in the laminar region</li>
<li>0 &le; Re &le; 1e6, 0.6 &le; Pr &le; 100, d/L &le; 1</li>
<li>The correlation holds for non-circular pipes only in the turbulent region. Use d_h=4*A/P as characteristic length.</li>
</ul>
The correlation takes into account the spatial position along the pipe flow, which changes discontinuously at flow reversal. However, the heat transfer coefficient itself is continuous around zero flow rate, but not its derivative.
<h4><font color=\"#008000\">References</font></h4>
 
<dl><dt>Verein Deutscher Ingenieure (1997):</dt>
    <dd><b>VDI W&auml;rmeatlas</b>.
         Springer Verlag, Ed. 8, 1997.</dd>
</dl>
</html>"));
    end PipeHT_LamTurb_local;
  end HeatTransfer;
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Pipes;
