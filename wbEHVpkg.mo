package wbEHVpkg "Package containing basic EV models"
  package EV
    model FirstEV "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = 0.26, S = 2.2, fc = 0.014, m = mass.m, rho(displayUnit = "kg/m3") = 1.226) annotation(
        Placement(visible = true, transformation(origin = {82, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 4) annotation(
        Placement(visible = true, transformation(origin = {-14, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
        Placement(visible = true, transformation(extent = {{-82, 10}, {-62, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "Sort1.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 1000, yMax = 100000.0) annotation(
        Placement(visible = true, transformation(extent = {{-118, 14}, {-98, 34}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{34, 10}, {54, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(extent = {{4, 10}, {24, 30}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 1.5) annotation(
        Placement(visible = true, transformation(extent = {{-54, 10}, {-34, 30}}, rotation = 0)));
    equation
      connect(torque.tau, driver.tauRef) annotation(
        Line(points = {{-84, 20}, {-90.5, 20}, {-90.5, 24}, {-97, 24}}, color = {0, 0, 127}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-108, 13}, {-108, -36}, {60, -36}, {60, -21}}, color = {0, 0, 127}));
      connect(mass.flange_b, dragF.flange) annotation(
        Line(points = {{54, 20}, {82, 20}, {82, -16}}, color = {0, 127, 0}));
      connect(velSens.flange, mass.flange_b) annotation(
        Line(points = {{60, 0}, {60, 20}, {54, 20}}, color = {0, 127, 0}));
      connect(mass.flange_a, wheel.flangeT) annotation(
        Line(points = {{34, 20}, {24, 20}}, color = {0, 127, 0}));
      connect(inertia.flange_a, torque.flange) annotation(
        Line(points = {{-54, 20}, {-58, 20}, {-62, 20}}));
      connect(inertia.flange_b, gear.flange_a) annotation(
        Line(points = {{-34, 20}, {-24, 20}}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-4, 20}, {-4, 20}, {4, 20}}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Very basic introductory EV model</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -40}, {100, 40}}, preserveAspectRatio = false), graphics = {Rectangle(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-84, 36}, {-24, 4}}), Text(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-82, 2}, {-26, -4}}, textString = "electric drive")}),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 200, Tolerance = 0.0001, Interval = 0.1),
        __OpenModelica_simulationFlags(jacobian = "", s = "dassl", lv = "LOG_STATS"),
        __OpenModelica_commandLineOptions = "");
    end FirstEV;

    model EVdata "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragF(m = mass.m, rho(displayUnit = "kg/m3") = data.rho, S = data.S, fc = data.fc, Cx = data.Cx) annotation(
        Placement(visible = true, transformation(origin = {82, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = data.ratio) annotation(
        Placement(visible = true, transformation(origin = {-14, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
        Placement(visible = true, transformation(extent = {{-82, 10}, {-62, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, yMax = 100000.0, k = data.kContr) annotation(
        Placement(visible = true, transformation(extent = {{-118, 10}, {-98, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = data.m) annotation(
        Placement(visible = true, transformation(extent = {{34, 10}, {54, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = data.radius) annotation(
        Placement(visible = true, transformation(extent = {{4, 10}, {24, 30}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J = data.J) annotation(
        Placement(visible = true, transformation(extent = {{-54, 10}, {-34, 30}}, rotation = 0)));
      VehicleData.Car data annotation(
        Placement(transformation(extent = {{30, 40}, {50, 60}})));
    equation
      connect(mass.flange_b, dragF.flange) annotation(
        Line(points = {{54, 20}, {82, 20}, {82, -16}}, color = {0, 127, 0}));
      connect(velSens.flange, mass.flange_b) annotation(
        Line(points = {{60, 0}, {60, 20}, {54, 20}}, color = {0, 127, 0}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-108, 8.8}, {-108, -36}, {60, -36}, {60, -21}}, color = {0, 0, 127}));
      connect(mass.flange_a, wheel.flangeT) annotation(
        Line(points = {{34, 20}, {24, 20}}, color = {0, 127, 0}));
      connect(inertia.flange_a, torque.flange) annotation(
        Line(points = {{-54, 20}, {-58, 20}, {-62, 20}}));
      connect(inertia.flange_b, gear.flange_a) annotation(
        Line(points = {{-34, 20}, {-24, 20}}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-4, 20}, {-4, 20}, {4, 20}}));
      connect(torque.tau, driver.tauRef) annotation(
        Line(points = {{-84, 20}, {-97, 20}}, color = {0, 0, 127}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Very basic introductory EV model</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -40}, {100, 60}}, preserveAspectRatio = false), graphics = {Rectangle(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-84, 36}, {-24, 4}}), Text(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-82, 2}, {-26, -4}}, textString = "electric drive")}),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 200, Tolerance = 0.0001, Interval = 0.1),
        __OpenModelica_simulationFlags(jacobian = "", s = "dassl", lv = "LOG_STATS"),
        __OpenModelica_commandLineOptions = "");
    end EVdata;

    model FirstEVpow "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = 0.26, S = 2.2, fc = 0.014, m = mass.m, rho (displayUnit = "kg/m3") = 1.226) annotation(
        Placement(visible = true, transformation(origin = {100, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 4) annotation(
        Placement(visible = true, transformation(origin = {-8, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
        Placement(visible = true, transformation(extent = {{-76, 10}, {-56, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "Sort1.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 1000, yMax = 100000.0) annotation(
        Placement(visible = true, transformation(extent = {{-112, 8}, {-92, 28}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{54, 10}, {74, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {78, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(extent = {{10, 10}, {30, 30}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 1.5) annotation(
        Placement(visible = true, transformation(extent = {{-48, 10}, {-28, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP2 annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {100, 6})));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP1 annotation(
        Placement(transformation(extent = {{36, 12}, {52, 28}})));
    equation
      connect(torque.tau, driver.tauRef) annotation(
        Line(points = {{-78, 20}, {-91, 20}}, color = {0, 0, 127}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-102, 8.8}, {-102, -36}, {78, -36}, {78, -21}}, color = {0, 0, 127}));
      connect(inertia.flange_b, gear.flange_a) annotation(
        Line(points = {{-28, 20}, {-18, 20}}));
      connect(inertia.flange_a, torque.flange) annotation(
        Line(points = {{-48, 20}, {-52, 20}, {-56, 20}}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{2, 20}, {2, 20}, {10, 20}}));
      connect(velSens.flange, mass.flange_b) annotation(
        Line(points = {{78, 4.44089e-016}, {78, 20}, {74, 20}}, color = {0, 127, 0}));
      connect(pP2.flange_b, dragF.flange) annotation(
        Line(points = {{100, -4}, {100, -12}}, color = {0, 127, 0}));
      connect(pP2.flange_a, mass.flange_b) annotation(
        Line(points = {{100, 16}, {100, 20}, {74, 20}}, color = {0, 127, 0}));
      connect(mass.flange_a, pP1.flange_b) annotation(
        Line(points = {{54, 20}, {54, 20}, {52, 20}}, color = {0, 127, 0}));
      connect(wheel.flangeT, pP1.flange_a) annotation(
        Line(points = {{30, 20}, {34, 20}, {36, 20}}, color = {0, 127, 0}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Very basic introductory EV model with power measurements.</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -40}, {120, 40}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(extent = {{-84, 36}, {-24, 4}}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash), Text(extent = {{-82, 2}, {-26, -4}}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, textString = "electric drive")}),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 200, Tolerance = 0.0001, Interval = 0.1),
        __OpenModelica_simulationFlags(jacobian = "", s = "dassl", lv = "LOG_STATS"));
    end FirstEVpow;

    model MBEVdata "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      Modelica.SIunits.Energy enBatDel, enDTdel, enP1del, enBattLoss, enBraking;
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = data.ratio) annotation(
        Placement(visible = true, transformation(origin = {-20, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = data.kContr, yMax = 100000.0) annotation(
        Placement(visible = true, transformation(extent = {{-116, -10}, {-96, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.5715) annotation(
        Placement(visible = true, transformation(extent = {{-4, 4}, {16, 24}}, rotation = 0)));
      wbEHPTlib.MapBased.OneFlange eleDrive(J = data.J, effTableName = "effTable", mapsFileName = "EVmaps.txt", mapsOnFile = true, powMax = 22e3, tauMax = 200, wMax = 1000) "Electric Drive" annotation(
        Placement(visible = true, transformation(extent = {{-74, 6}, {-54, 24}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1 batt1(SOCInit = 0.7, QCellNom = 100 * 3600, ns = 100) annotation(
        Placement(transformation(extent = {{-112, 34}, {-92, 54}})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(extent = {{-84, -20}, {-64, 0}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP1 annotation(
        Placement(visible = true, transformation(origin = {32, 14}, extent = {{-6, -8}, {6, 8}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = data.m) annotation(
        Placement(visible = true, transformation(extent = {{56, 4}, {76, 24}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP2 annotation(
        Placement(visible = true, transformation(origin = {98, 4}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {68, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = data.Cx, S = data.S, fc = data.fc, m = data.m, rho = data.rho) annotation(
        Placement(visible = true, transformation(origin = {98, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Mechanics.Translational.Sources.Force brake annotation(
        Placement(visible = true, transformation(origin = {32, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sensors.TorqueSensor edTau annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 0, origin = {-42, 14})));
      Modelica.Blocks.Math.Add add(k1 = -1) annotation(
        Placement(transformation(extent = {{-42, -26}, {-30, -14}})));
      Modelica.Blocks.Math.Gain tqToForce(k = gear.ratio * wheel.radius) annotation(
        Placement(visible = true, transformation(extent = {{0, -26}, {12, -14}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter cutNeg(limitsAtInit = true, uMax = 0, uMin = -Modelica.Constants.inf) annotation(
        Placement(visible = true, transformation(origin = {-14, -20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  VehicleData.Car data annotation(
        Placement(visible = true, transformation(origin = {70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(batt1.n, eleDrive.pin_n) annotation(
        Line(points = {{-91.9, 38}, {-80, 38}, {-80, 10}, {-74, 10}}, color = {0, 0, 255}));
      connect(eleDrive.pin_n, ground.p) annotation(
        Line(points = {{-74, 10}, {-74, 10}, {-74, 0}, {-74, 0}}, color = {0, 0, 255}));
      connect(eleDrive.tauRef, driver.tauRef) annotation(
        Line(points = {{-75.4, 14}, {-86, 14}, {-86, 0}, {-95, 0}}, color = {0, 0, 127}));
      connect(batt1.p, eleDrive.pin_p) annotation(
        Line(points = {{-92, 50}, {-74, 50}, {-74, 18}}, color = {0, 0, 255}));
      connect(edTau.flange_a, eleDrive.flange_a) annotation(
        Line(points = {{-50, 14}, {-54, 14}}));
      connect(tqToForce.y, brake.f) annotation(
        Line(points = {{12.6, -20}, {18, -20}, {18, -20}, {20, -20}}, color = {0, 0, 127}));
      connect(cutNeg.y, tqToForce.u) annotation(
        Line(points = {{-7.4, -20}, {-2, -20}, {-2, -20}, {-1.2, -20}}, color = {0, 0, 127}));
      connect(cutNeg.u, add.y) annotation(
        Line(points = {{-21.2, -20}, {-30, -20}, {-30, -20}, {-29.4, -20}}, color = {0, 0, 127}));
      connect(brake.flange, pP1.flange_b) annotation(
        Line(points = {{42, -20}, {46, -20}, {46, 14}, {38, 14}}, color = {0, 127, 0}));
      connect(velSens.flange, pP2.flange_a) annotation(
        Line(points = {{78, -42}, {78, 14}, {98, 14}}, color = {0, 127, 0}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-106, -11.2}, {-106, -42}, {57, -42}}, color = {0, 0, 127}));
      connect(mass.flange_a, pP1.flange_b) annotation(
        Line(points = {{56, 14}, {38, 14}}, color = {0, 127, 0}));
      connect(pP2.flange_a, mass.flange_b) annotation(
        Line(points = {{98, 14}, {76, 14}}, color = {0, 127, 0}));
      connect(pP1.flange_a, wheel.flangeT) annotation(
        Line(points = {{26, 14}, {16, 14}}, color = {0, 127, 0}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-10, 14}, {-4, 14}}));
      connect(dragF.flange, pP2.flange_b) annotation(
        Line(points = {{98, -14}, {98, -6}}, color = {0, 127, 0}));
      der(enBatDel) = (batt1.p.v - batt1.n.v) * batt1.n.i;
      der(enDTdel) = eleDrive.powSensor.power;
      der(enP1del) = pP1.power;
      der(enBattLoss) = batt1.powerLoss;
      der(enBraking) = if pP1.power > 0 then 0 else -pP1.power;
      connect(add.u2, driver.tauRef) annotation(
        Line(points = {{-43.2, -23.6}, {-86, -23.6}, {-86, 0}, {-95, 0}}, color = {0, 0, 127}));
      connect(edTau.flange_b, gear.flange_a) annotation(
        Line(points = {{-34, 14}, {-30, 14}, {-30, 14}}, color = {0, 0, 0}));
      connect(add.u1, edTau.tau) annotation(
        Line(points = {{-43.2, -16.4}, {-48.4, -16.4}, {-48.4, 5.2}}, color = {0, 0, 127}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Simple map-based EV model with battery.</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -60}, {120, 60}}, preserveAspectRatio = false)),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 0.1),
        __OpenModelica_commandLineOptions = "");
    end MBEVdata;


    model FirstEVAngle "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 6) annotation(
        Placement(visible = true, transformation(origin = {-14, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation(
        Placement(visible = true, transformation(extent = {{-82, 10}, {-62, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 1000, yMax = 100000.0, CycleFileName = "TestAngle.txt") annotation(
        Placement(visible = true, transformation(extent = {{-118, 10}, {-98, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{34, 10}, {54, 30}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.5715) annotation(
        Placement(visible = true, transformation(extent = {{4, 10}, {24, 30}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 5) annotation(
        Placement(visible = true, transformation(extent = {{-54, 10}, {-34, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForceAngle dragF1(Cx = 0.65, S = 6.0, fc = 0.013, m = mass.m, rho = 1.226, DataFileName = "Angle1.txt") annotation(
        Placement(visible = true, transformation(origin = {82, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(velSens.flange, mass.flange_b) annotation(
        Line(points = {{60, 0}, {60, 20}, {54, 20}}, color = {0, 127, 0}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-108, 8.8}, {-108, -36}, {60, -36}, {60, -21}}, color = {0, 0, 127}));
      connect(mass.flange_a, wheel.flangeT) annotation(
        Line(points = {{34, 20}, {24, 20}}, color = {0, 127, 0}));
      connect(inertia.flange_a, torque.flange) annotation(
        Line(points = {{-54, 20}, {-58, 20}, {-62, 20}}));
      connect(inertia.flange_b, gear.flange_a) annotation(
        Line(points = {{-34, 20}, {-24, 20}}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-4, 20}, {-4, 20}, {4, 20}}));
      connect(torque.tau, driver.tauRef) annotation(
        Line(points = {{-84, 20}, {-97, 20}}, color = {0, 0, 127}));
      connect(dragF1.flange, mass.flange_b) annotation(
        Line(points = {{82, -14}, {82, 20}, {54, 20}}, color = {0, 127, 0}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Very basic introductory EV model</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -40}, {100, 40}}, preserveAspectRatio = false), graphics = {Rectangle(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-84, 36}, {-24, 4}}), Text(origin = {-6, 0}, lineColor = {28, 108, 200}, pattern = LinePattern.Dash, extent = {{-82, 2}, {-26, -4}}, textString = "electric drive")}),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StopTime = 200, Interval = 0.1),
        __OpenModelica_simulationFlags(jacobian = "", s = "dassl", lv = "LOG_STATS"),
        __OpenModelica_commandLineOptions = "");
    end FirstEVAngle;

    model MBEV "Simulates a very basic Electric Vehicle"
      import Modelica;
      extends Modelica.Icons.Example;
      Modelica.SIunits.Energy enBatDel, enDTdel, enP1del, enBattLoss, enBraking;
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 6) annotation(
        Placement(visible = true, transformation(origin = {-20, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, yMax = 100000.0, k = 100) annotation(
        Placement(visible = true, transformation(extent = {{-116, -10}, {-96, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.5715) annotation(
        Placement(visible = true, transformation(extent = {{-4, 4}, {16, 24}}, rotation = 0)));
      wbEHPTlib.MapBased.OneFlange eleDrive(J = 0.25, effTableName = "effTable", mapsFileName = "EVmaps.txt", tauMax = 150, wMax = 500) "Electric Drive" annotation(
        Placement(visible = true, transformation(extent = {{-74, 6}, {-54, 24}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1 batt1(SOCInit = 0.7, QCellNom = 100 * 3600, ns = 100) annotation(
        Placement(transformation(extent = {{-112, 34}, {-92, 54}})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(extent = {{-84, -20}, {-64, 0}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP1 annotation(
        Placement(visible = true, transformation(origin = {32, 14}, extent = {{-6, -8}, {6, 8}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{56, 4}, {76, 24}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor pP2 annotation(
        Placement(visible = true, transformation(origin = {98, 4}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor velSens annotation(
        Placement(visible = true, transformation(origin = {68, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragF(Cx = 0.26, rho = 1.226, S = 2.2, fc = 0.014, m = mass.m) annotation(
        Placement(visible = true, transformation(origin = {98, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Mechanics.Translational.Sources.Force brake annotation(
        Placement(visible = true, transformation(origin = {32, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sensors.TorqueSensor edTau annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 0, origin = {-42, 14})));
      Modelica.Blocks.Math.Add add(k1 = -1) annotation(
        Placement(transformation(extent = {{-42, -26}, {-30, -14}})));
      Modelica.Blocks.Math.Gain tqToForce(k = gear.ratio * wheel.radius) annotation(
        Placement(visible = true, transformation(extent = {{0, -26}, {12, -14}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter cutNeg(limitsAtInit = true, uMax = 0, uMin = -Modelica.Constants.inf) annotation(
        Placement(visible = true, transformation(origin = {-14, -20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    equation
      connect(batt1.n, eleDrive.pin_n) annotation(
        Line(points = {{-91.9, 38}, {-80, 38}, {-80, 10}, {-74, 10}}, color = {0, 0, 255}));
      connect(eleDrive.pin_n, ground.p) annotation(
        Line(points = {{-74, 10}, {-74, 10}, {-74, 0}, {-74, 0}}, color = {0, 0, 255}));
      connect(eleDrive.tauRef, driver.tauRef) annotation(
        Line(points = {{-75.4, 14}, {-86, 14}, {-86, 0}, {-95, 0}}, color = {0, 0, 127}));
      connect(batt1.p, eleDrive.pin_p) annotation(
        Line(points = {{-92, 50}, {-74, 50}, {-74, 18}}, color = {0, 0, 255}));
      connect(edTau.flange_a, eleDrive.flange_a) annotation(
        Line(points = {{-50, 14}, {-54, 14}}));
      connect(tqToForce.y, brake.f) annotation(
        Line(points = {{12.6, -20}, {18, -20}, {18, -20}, {20, -20}}, color = {0, 0, 127}));
      connect(cutNeg.y, tqToForce.u) annotation(
        Line(points = {{-7.4, -20}, {-2, -20}, {-2, -20}, {-1.2, -20}}, color = {0, 0, 127}));
      connect(cutNeg.u, add.y) annotation(
        Line(points = {{-21.2, -20}, {-30, -20}, {-30, -20}, {-29.4, -20}}, color = {0, 0, 127}));
      connect(brake.flange, pP1.flange_b) annotation(
        Line(points = {{42, -20}, {46, -20}, {46, 14}, {38, 14}}, color = {0, 127, 0}));
      connect(velSens.flange, pP2.flange_a) annotation(
        Line(points = {{78, -42}, {78, 14}, {98, 14}}, color = {0, 127, 0}));
      connect(driver.V, velSens.v) annotation(
        Line(points = {{-106, -11.2}, {-106, -42}, {57, -42}}, color = {0, 0, 127}));
      connect(mass.flange_a, pP1.flange_b) annotation(
        Line(points = {{56, 14}, {38, 14}}, color = {0, 127, 0}));
      connect(pP2.flange_a, mass.flange_b) annotation(
        Line(points = {{98, 14}, {76, 14}}, color = {0, 127, 0}));
      connect(pP1.flange_a, wheel.flangeT) annotation(
        Line(points = {{26, 14}, {16, 14}}, color = {0, 127, 0}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-10, 14}, {-4, 14}}));
      connect(dragF.flange, pP2.flange_b) annotation(
        Line(points = {{98, -14}, {98, -6}}, color = {0, 127, 0}));
      der(enBatDel) = (batt1.p.v - batt1.n.v) * batt1.n.i;
      der(enDTdel) = eleDrive.powSensor.power;
      der(enP1del) = pP1.power;
      der(enBattLoss) = batt1.powerLoss;
      der(enBraking) = if pP1.power > 0 then 0 else -pP1.power;
      connect(add.u2, driver.tauRef) annotation(
        Line(points = {{-43.2, -23.6}, {-86, -23.6}, {-86, 0}, {-95, 0}}, color = {0, 0, 127}));
      connect(edTau.flange_b, gear.flange_a) annotation(
        Line(points = {{-34, 14}, {-30, 14}, {-30, 14}}, color = {0, 0, 0}));
      connect(add.u1, edTau.tau) annotation(
        Line(points = {{-43.2, -16.4}, {-48.4, -16.4}, {-48.4, 5.2}}, color = {0, 0, 127}));
      annotation(
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>Simple map-based EV model with battery.</p>
</html>"),
        Commands,
        Diagram(coordinateSystem(extent = {{-120, -60}, {120, 60}}, preserveAspectRatio = false)),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 0.1),
        __OpenModelica_commandLineOptions = "");
    end MBEV;


  end EV;

  package sHEV
    model SHEVpowerFilt "Ice, Generator, DriveTrain, all map-based"
      //€
      extends Modelica.Icons.Example;
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 0, origin = {-42, -2})));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.473) annotation(
        Placement(transformation(extent = {{-52, -46}, {-38, -32}})));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 10) annotation(
        Placement(transformation(extent = {{-78, -46}, {-64, -32}})));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 14000) annotation(
        Placement(transformation(extent = {{-8, -48}, {10, -30}})));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powProp annotation(
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 0, origin = {-23, -39})));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powDrag annotation(
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 0, origin = {45, -39})));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(S = 6.5, fc = 0.01, Cx = 0.65, m = mass.m) annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 90, origin = {80, -48})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor1 annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 270, origin = {26, -54})));
      Modelica.Blocks.Continuous.FirstOrder powFilt(y_start = 20e3, T = 500) annotation(
        Placement(visible = true, transformation(extent = {{12, 64}, {-4, 80}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiter(uMax = 100e3, uMin = 0) annotation(
        Placement(visible = true, transformation(extent = {{-14, 64}, {-30, 80}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "Sort1.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 200.0, yMax = 2e3) annotation(
        Placement(visible = true, transformation(extent = {{-98, 76}, {-78, 96}}, rotation = 0)));
      wbEHPTlib.MapBased.Genset genset(OptiTime = 1, mapsFileName = "SHEVmaps.txt", maxGenW = 300, maxTau = 400, maxPow = 45000) annotation(
        Placement(transformation(extent = {{-80, 8}, {-50, 38}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1 battery(QCellNom = 200 * 3600, ICellMax = 1000, R0Cell = 0.35E-3, efficiency = 0.9, iCellEfficiency = 200, ns = 200) annotation(
        Placement(transformation(extent = {{0, 18}, {20, 38}})));
      wbEHPTlib.MapBased.OneFlange drive( effTableName = "motEffTable", mapsFileName = "SHEVmaps.txt", mapsOnFile = true, powMax = 150e3,tauMax = 1000, wMax = 3000) annotation(
        Placement(visible = true, transformation(extent = {{40, 42}, {60, 22}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor drivePow annotation(
        Placement(transformation(extent = {{66, 40}, {86, 20}})));
      Modelica.Electrical.Analog.Sensors.PowerSensor gsPow annotation(
        Placement(transformation(extent = {{-32, 22}, {-12, 42}})));
    equation
      connect(drive.pin_n, battery.p) annotation(
        Line(points = {{40, 37.5556}, {36.1, 37.5556}, {36.1, 42}, {20, 42}, {20, 34}}, color = {0, 0, 255}));
      connect(gsPow.nv, drive.pin_p) annotation(
        Line(points = {{-22, 22}, {-22, 14}, {40, 14}, {40, 28.6667}}, color = {0, 0, 255}));
      connect(drivePow.flange_a, drive.flange_a) annotation(
        Line(points = {{66, 30}, {61, 30}, {61, 33.1111}, {60, 33.1111}}));
      connect(battery.n, drive.pin_p) annotation(
        Line(points = {{20.1, 22}, {24, 22}, {24, 14}, {40, 14}, {40, 28.6667}}, color = {0, 0, 255}));
      connect(genset.pin_n, drive.pin_p) annotation(
        Line(points = {{-49.7, 14}, {40, 14}, {40, 28.6667}}, color = {0, 0, 255}));
      connect(drive.tauRef, driver.tauRef) annotation(
        Line(points = {{38.6, 33.1111}, {38.6, 34}, {32, 34}, {32, 86}, {-77, 86}}, color = {0, 0, 127}));
      connect(drivePow.power, powFilt.u) annotation(
        Line(points = {{68, 41}, {68, 72}, {13.6, 72}}, color = {0, 0, 127}));
      connect(powFilt.y, limiter.u) annotation(
        Line(points = {{-4.8, 72}, {-12.4, 72}}, color = {0, 0, 127}));
      connect(limiter.y, genset.powRef) annotation(
        Line(points = {{-30.8, 72}, {-55.85, 72}, {-55.85, 40.25}}, color = {0, 0, 127}));
      connect(speedSensor1.v, driver.V) annotation(
        Line(points = {{26, -62.8}, {26, -74}, {-94, -74}, {-94, 68}, {-88, 68}, {-88, 74.8}}, color = {0, 0, 127}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-64, -39}, {-52, -39}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(mass.flange_a, powProp.flange_b) annotation(
        Line(points = {{-8, -39}, {-16, -39}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powProp.flange_a, wheel.flangeT) annotation(
        Line(points = {{-30, -39}, {-38, -39}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powDrag.flange_a, mass.flange_b) annotation(
        Line(points = {{38, -39}, {10, -39}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powDrag.flange_b, dragForce.flange) annotation(
        Line(points = {{52, -39}, {66, -39}, {66, -40}, {80, -40}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(speedSensor1.flange, mass.flange_b) annotation(
        Line(points = {{26, -46}, {26, -39}, {10, -39}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(genset.pin_n, ground1.p) annotation(
        Line(points = {{-49.7, 14}, {-42, 14}, {-42, 6}}, color = {0, 0, 255}, smooth = Smooth.None));
      connect(drivePow.flange_b, gear.flange_a) annotation(
        Line(points = {{86, 30}, {90, 30}, {90, -22}, {-86, -22}, {-86, -39}, {-78, -39}}, color = {0, 0, 0}));
      connect(gsPow.pc, genset.pin_p) annotation(
        Line(points = {{-32, 32}, {-36, 32}, {-50, 32}}, color = {0, 0, 255}));
      connect(gsPow.nc, battery.p) annotation(
        Line(points = {{-12, 32}, {-8, 32}, {-8, 48}, {28, 48}, {28, 34}, {20, 34}}, color = {0, 0, 255}));
      connect(gsPow.pv, gsPow.pc) annotation(
        Line(points = {{-22, 42}, {-32, 42}, {-32, 32}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -80}, {100, 100}}), graphics = {Rectangle(extent = {{-90, -28}, {94, -70}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Rectangle(extent = {{-90, 52}, {96, -10}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Rectangle(extent = {{-60, 96}, {94, 58}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Text(extent = {{68, 74}, {94, 66}}, lineColor = {255, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "EMS"), Text(extent = {{-96, -60}, {-44, -68}}, lineColor = {255, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "MechProp"), Text(extent = {{12, 0}, {58, -8}}, lineColor = {255, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "PowerTrain")}),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1500, StartTime = 0, Tolerance = 1e-06, Interval = 0.75),
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>This is a SHEV model which has an Energy Management System able to control the power flow with basic logic: requests the ICE to deliver the average load power.</p>
</html>"));
    end SHEVpowerFilt;

    model SHEVpowerFiltSoc "Ice, Generator, DriveTrain, all map-based"
      //€
      extends Modelica.Icons.Example;
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(transformation(extent = {{-8, 8}, {8, -8}}, rotation = 0, origin = {-42, -4})));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.473) annotation(
        Placement(transformation(extent = {{-52, -50}, {-40, -38}})));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 10) annotation(
        Placement(transformation(extent = {{-78, -50}, {-66, -38}})));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 14000) annotation(
        Placement(transformation(extent = {{-8, -52}, {10, -34}})));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powProp annotation(
        Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 0, origin = {-24, -44})));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powDrag annotation(
        Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 0, origin = {45, -43})));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(S = 6.5, fc = 0.01, Cx = 0.65, m = mass.m) annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 90, origin = {80, -52})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor1 annotation(
        Placement(transformation(extent = {{-8, -8}, {8, 8}}, rotation = 270, origin = {26, -58})));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 200.0, yMax = 2e3) annotation(
        Placement(visible = true, transformation(extent = {{-100, 94}, {-80, 114}}, rotation = 0)));
      wbEHPTlib.MapBased.Genset genset(OptiTime = 1, mapsFileName = "SHEVmaps.txt") annotation(
        Placement(transformation(extent = {{-80, -18}, {-50, 12}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1 battery(ICellMax = 1000, QCellNom = 200 * 3600, R0Cell = 0.35E-3, SOCInit = 0.6, efficiency = 0.9, iCellEfficiency = 200, ns = 200) annotation(
        Placement(transformation(extent = {{0, -8}, {20, 12}})));
      wbEHPTlib.MapBased.OneFlange gen(tauMax = 1000, powMax = 150e3, wMax = 3000, mapsFileName = "SHEVmaps.txt", effTableName = "motEffTable") annotation(
        Placement(visible = true, transformation(extent = {{40, 16}, {60, -4}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor drivePow annotation(
        Placement(transformation(extent = {{66, 18}, {86, -2}})));
      Modelica.Electrical.Analog.Sensors.PowerSensor gsPow annotation(
        Placement(transformation(extent = {{-32, -4}, {-12, 16}})));
      Modelica.Blocks.Math.Feedback fbSOC annotation(
        Placement(transformation(extent = {{20, 40}, {0, 60}})));
      Modelica.Blocks.Sources.Constant socRef_(k = 0.6) annotation(
        Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = 180, origin = {38, 50})));
      Modelica.Blocks.Math.Gain socErrToPow(k = 3e6) annotation(
        Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 0, origin = {-24, 50})));
      Modelica.Blocks.Math.Add add annotation(
        Placement(visible = true, transformation(origin = {-6, 84}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.FirstOrder powFilt(y_start = 20e3, T = 500) annotation(
        Placement(visible = true, transformation(extent = {{34, 82}, {18, 98}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter toPowRef(uMax = 100e3, uMin = 0) annotation(
        Placement(visible = true, transformation(extent = {{-30, 76}, {-46, 92}}, rotation = 0)));
    equation
      connect(gsPow.nv, gen.pin_p) annotation(
        Line(points = {{-22, -4}, {-22, -12}, {32, -12}, {32, -1}, {40, -1}, {40, 2.66667}}, color = {0, 0, 255}));
      connect(battery.n, gen.pin_p) annotation(
        Line(points = {{20.1, -4}, {24, -4}, {24, -12}, {32, -12}, {32, -1}, {40, -1}, {40, 2.66667}}, color = {0, 0, 255}));
      connect(genset.pin_n, gen.pin_p) annotation(
        Line(points = {{-49.7, -12}, {32, -12}, {32, -1}, {40, -1}, {40, 2.66667}}, color = {0, 0, 255}));
      connect(gen.tauRef, driver.tauRef) annotation(
        Line(points = {{38.6, 7.11111}, {38.6, 6}, {34, 6}, {34, 28}, {62, 28}, {62, 104}, {-79, 104}}, color = {0, 0, 127}));
      connect(gen.pin_n, battery.p) annotation(
        Line(points = {{40, 11.5556}, {28, 11.5556}, {28, 8}, {20, 8}}, color = {0, 0, 255}));
      connect(drivePow.flange_a, gen.flange_a) annotation(
        Line(points = {{66, 8}, {61, 8}, {61, 7.11111}, {60, 7.11111}}));
      connect(socErrToPow.y, add.u2) annotation(
        Line(points = {{-35, 50}, {-46, 50}, {-46, 68}, {12, 68}, {12, 78}, {6, 78}}, color = {0, 0, 127}));
      connect(add.y, toPowRef.u) annotation(
        Line(points = {{-17, 84}, {-28.4, 84}}, color = {0, 0, 127}));
      connect(toPowRef.y, genset.powRef) annotation(
        Line(points = {{-46.8, 84}, {-56, 84}, {-56, 14.25}, {-55.85, 14.25}}, color = {0, 0, 127}));
      connect(powFilt.y, add.u1) annotation(
        Line(points = {{17.2, 90}, {6, 90}}, color = {0, 0, 127}));
      connect(drivePow.power, powFilt.u) annotation(
        Line(points = {{68, 19}, {68, 90}, {35.6, 90}}, color = {0, 0, 127}));
      connect(speedSensor1.v, driver.V) annotation(
        Line(points = {{26, -66.8}, {26, -78}, {-94, -78}, {-94, 44}, {-90, 44}, {-90, 92.8}}, color = {0, 0, 127}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-66, -44}, {-52, -44}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(mass.flange_a, powProp.flange_b) annotation(
        Line(points = {{-8, -43}, {-12, -43}, {-12, -44}, {-18, -44}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powProp.flange_a, wheel.flangeT) annotation(
        Line(points = {{-30, -44}, {-40, -44}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powDrag.flange_a, mass.flange_b) annotation(
        Line(points = {{38, -43}, {10, -43}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(powDrag.flange_b, dragForce.flange) annotation(
        Line(points = {{52, -43}, {66, -43}, {66, -44}, {80, -44}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(speedSensor1.flange, mass.flange_b) annotation(
        Line(points = {{26, -50}, {26, -43}, {10, -43}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(genset.pin_n, ground1.p) annotation(
        Line(points = {{-49.7, -12}, {-42, -12}}, color = {0, 0, 255}, smooth = Smooth.None));
      connect(drivePow.flange_b, gear.flange_a) annotation(
        Line(points = {{86, 8}, {90, 8}, {90, -26}, {-86, -26}, {-86, -44}, {-78, -44}}, color = {0, 0, 0}));
      connect(gsPow.pc, genset.pin_p) annotation(
        Line(points = {{-32, 6}, {-36, 6}, {-50, 6}}, color = {0, 0, 255}));
      connect(gsPow.nc, battery.p) annotation(
        Line(points = {{-12, 6}, {-8, 6}, {-8, 22}, {28, 22}, {28, 8}, {20, 8}}, color = {0, 0, 255}));
      connect(gsPow.pv, gsPow.pc) annotation(
        Line(points = {{-22, 16}, {-32, 16}, {-32, 6}}, color = {0, 0, 255}));
      connect(battery.SOC, fbSOC.u2) annotation(
        Line(points = {{-1, 2}, {-4, 2}, {-4, 34}, {-4, 42}, {10, 42}}, color = {0, 0, 127}));
      connect(fbSOC.u1, socRef_.y) annotation(
        Line(points = {{18, 50}, {27, 50}}, color = {0, 0, 127}));
      connect(socErrToPow.u, fbSOC.y) annotation(
        Line(points = {{-12, 50}, {1, 50}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -80}, {100, 120}}), graphics = {Rectangle(extent = {{-90, -32}, {94, -74}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Rectangle(extent = {{-90, 26}, {96, -22}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Rectangle(extent = {{-70, 114}, {96, 32}}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash), Text(extent = {{68, 92}, {94, 84}}, lineColor = {255, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "EMS"), Text(extent = {{-96, -64}, {-44, -72}}, lineColor = {255, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid, textString = "MechProp")}),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 2000, StartTime = 0, Tolerance = 1e-06, Interval = 4),
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>This is a SHEV model which has an Energy Management System able to control the power flow with:</p>
<p>- basic logic: requests the ICE to deliver the average load power </p>
<p>- additional logic: SOC loop to avoid SOC drift.</p>
</html>"));
    end SHEVpowerFiltSoc;

    model SHEV_OO "Ice, Generator, DriveTrain, all map-based"
      //€
      extends Modelica.Icons.Example;
      Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
        Placement(visible = true, transformation(origin = {22, -12}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.473) annotation(
        Placement(visible = true, transformation(origin = {-48, -42}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 10) annotation(
        Placement(visible = true, transformation(extent = {{-82, -52}, {-62, -32}}, rotation = 0)));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 14000) annotation(
        Placement(visible = true, transformation(extent = {{-2, -50}, {16, -32}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powProp annotation(
        Placement(visible = true, transformation(origin = {-21, -41}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PowerSensor powDrag annotation(
        Placement(visible = true, transformation(origin = {51, -41}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(S = 6.5, fc = 0.01, Cx = 0.65, m = mass.m) annotation(
        Placement(visible = true, transformation(origin = {76, -50}, extent = {{-8, -8}, {8, 8}}, rotation = 90)));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor speedSensor1 annotation(
        Placement(visible = true, transformation(origin = {26, -56}, extent = {{-8, -8}, {8, 8}}, rotation = 270)));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1 battery(ICellMax = 1000, QCellNom = 60 * 3600, R0Cell = 0.35E-3, SOCInit = 0.6, efficiency = 0.9, iCellEfficiency = 200, ns = 200) annotation(
        Placement(visible = true, transformation(origin = {-2, 26}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      wbEHPTlib.MapBased.OneFlange drive(tauMax = 1000, powMax = 150e3, wMax = 3000, mapsFileName = "SHEVmaps.txt", effTableName = "motEffTable") annotation(
        Placement(visible = true, transformation(extent = {{36, 16}, {56, -4}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor drivePow annotation(
        Placement(visible = true, transformation(extent = {{62, 14}, {82, -6}}, rotation = 0)));
      Modelica.Electrical.Analog.Sensors.PowerSensor gsPow annotation(
        Placement(visible = true, transformation(extent = {{-36, 0}, {-16, 20}}, rotation = 0)));
      wbEHPTlib.MapBased.ECUs.EMS ems(powHigh = 60e3, powLow = 30e3, powMax = 200e3, powPerSoc = 300e3, socRef = 0.75) annotation(
        Placement(visible = true, transformation(origin = {-20, 58}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      wbEHPTlib.MapBased.GensetOO genset(iceTauMaxReq = drive.tauMax, maxGensetTau = 700, maxIceTau = [80, 380; 100, 620; 120, 800; 170, 800; 220, 670; 230, 650; 240, 570], specConsumption = [0.0, 0, 80, 100, 120, 170, 220, 230, 240; 0, 999, 999, 999, 999, 999, 999, 999, 999; 100, 280, 280, 280, 280, 280, 280, 280, 280; 200, 215, 213, 219, 225, 231, 252, 255, 260; 300, 210, 208, 208, 207, 214, 228, 230, 237; 400, 213, 204, 204, 199, 204, 218, 218, 222; 500, 213, 205, 205, 198, 202, 213, 214, 215; 600, 213, 213, 213, 198, 205, 212, 213, 214; 700, 213, 213, 213, 198, 204, 210, 211, 214; 800, 213, 213, 213, 213, 213, 210, 211, 214], optiTable = [0, 800; 20000, 850; 40000, 1100; 60000, 1250; 80000, 1280; 100000, 1340; 120000, 1400; 140000, 1650; 160000, 2130], maxPow = 45e3, mapsFileName = "SHEVmaps.txt") annotation(
        Placement(visible = true, transformation(extent = {{-84, -14}, {-54, 16}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "Sort1.txt", extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, k = 200.0, yMax = 2e3) annotation(
        Placement(visible = true, transformation(origin = {86, 60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(battery.n, drive.pin_p) annotation(
        Line(points = {{-8, 15.9}, {-8, -8}, {10, -8}, {10, 2.66667}, {36, 2.66667}}, color = {0, 0, 255}));
      connect(gsPow.nv, drive.pin_p) annotation(
        Line(points = {{-26, 0}, {-26, -8}, {10, -8}, {10, 2.66667}, {36, 2.66667}}, color = {0, 0, 255}));
      connect(genset.pin_n, drive.pin_p) annotation(
        Line(points = {{-53.7, -8}, {-32, -8}, {10, -8}, {10, 2.66667}, {36, 2.66667}}, color = {0, 0, 255}));
      connect(ground1.p, drive.pin_p) annotation(
        Line(points = {{22, -4}, {22, 2.66667}, {36, 2.66667}}, color = {0, 0, 255}));
      connect(drivePow.flange_a, drive.flange_a) annotation(
        Line(points = {{62, 4}, {58, 4}, {58, 7.11111}, {56, 7.11111}}));
      connect(drive.pin_n, gsPow.nc) annotation(
        Line(points = {{36, 11.5556}, {22, 11.5556}, {22, 10}, {-16, 10}}, color = {0, 0, 255}));
      connect(driver.tauRef, drive.tauRef) annotation(
        Line(points = {{75, 60}, {34.6, 60}, {34.6, 7.11111}}, color = {0, 0, 127}));
      connect(speedSensor1.v, driver.V) annotation(
        Line(points = {{26, -64.8}, {26, -74}, {98, -74}, {98, 44}, {86, 44}, {86, 48.8}}, color = {0, 0, 127}));
      connect(gsPow.pc, genset.pin_p) annotation(
        Line(points = {{-36, 10}, {-54, 10}}, color = {0, 0, 255}));
      connect(ems.pcPowReq, genset.powRef) annotation(
        Line(points = {{-30.8, 52}, {-60, 52}, {-60, 18.4}}, color = {0, 0, 127}));
      connect(ems.on, genset.ON) annotation(
        Line(points = {{-30.8, 64}, {-78, 64}, {-78, 18.4}}, color = {255, 0, 255}));
      connect(gear.flange_b, wheel.flangeR) annotation(
        Line(points = {{-62, -42}, {-56, -42}}));
      connect(powProp.flange_a, wheel.flangeT) annotation(
        Line(points = {{-30, -41}, {-35, -41}, {-35, -42}, {-40, -42}}, color = {0, 127, 0}));
      connect(mass.flange_a, powProp.flange_b) annotation(
        Line(points = {{-2, -41}, {-12, -41}}, color = {0, 127, 0}));
      connect(battery.SOC, ems.soc) annotation(
        Line(points = {{-2, 37}, {-2, 52}, {-8, 52}}, color = {0, 0, 127}));
      connect(ems.edPow, drivePow.power) annotation(
        Line(points = {{-8, 64}, {38, 64}, {38, 64}, {64, 64}, {64, 15}}, color = {0, 0, 127}));
      connect(gsPow.pv, gsPow.pc) annotation(
        Line(points = {{-26, 20}, {-36, 20}, {-36, 10}}, color = {0, 0, 255}));
      connect(powDrag.flange_b, dragForce.flange) annotation(
        Line(points = {{60, -41}, {60, -42}, {76, -42}}, color = {0, 127, 0}));
      connect(drivePow.flange_b, gear.flange_a) annotation(
        Line(points = {{82, 4}, {88, 4}, {88, -18}, {-86, -18}, {-86, -42}, {-82, -42}}));
      connect(speedSensor1.flange, mass.flange_b) annotation(
        Line(points = {{26, -48}, {26, -41}, {16, -41}}, color = {0, 127, 0}));
      connect(powDrag.flange_a, mass.flange_b) annotation(
        Line(points = {{42, -41}, {42, -41}, {16, -41}}, color = {0, 127, 0}));
      connect(battery.p, gsPow.nc) annotation(
        Line(points = {{4, 16}, {4, 10}, {-16, 10}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-100, -80}, {100, 80}}, initialScale = 0.1), graphics = {Rectangle(origin = {-2, 4}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash, extent = {{-90, -32}, {94, -74}}), Rectangle(origin = {-4, 0}, lineColor = {255, 0, 0}, pattern = LinePattern.Dash, extent = {{-88, 40}, {96, -22}}), Text(origin = {0, 4}, lineColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-96, -64}, {-44, -72}}, textString = "MechProp")}),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics),
        experiment(StopTime = 1800, __Dymola_NumberOfIntervals = 2000),
        experimentSetupOutput(derivatives = false),
        Documentation(info = "<html>
<p>This is a SHEV model which has an Energy Management System able to control the power flow with:</p>
<p>- basic logic: requests the ICE to deliver the average load power </p>
<p>- additional logic: SOC loop to avoid SOC drift</p>
<p>- further ON/OFF control to switch OFF the engine when the average power is too low to permit efficient operation.</p>
</html>"),
        __OpenModelica_commandLineOptions = "");
    end SHEV_OO;
  end sHEV;

  package PSD
    model PSecu1 "Full Power Split Device power train using Map-Based components"
      import Modelica.Constants.*;
      extends Modelica.Icons.Example;
      parameter Modelica.SIunits.AngularVelocity wIceStart = 50;
      Modelica.SIunits.Energy EbatDel "energy delivered by the battery";
      Modelica.SIunits.Energy EgenDelM "energy delivered by gen trough mechanical flange";
      Modelica.SIunits.Energy Eroad "mechanical energy absorbed by roas (friction & air)";
      Modelica.SIunits.Energy EiceDel "mechanical energy delivered by ice";
      Modelica.SIunits.Energy Emot, Emass;
      Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-50, 52})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
        Placement(transformation(extent = {{2, 42}, {22, 62}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor carVel annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {78, -12})));
      Modelica.Mechanics.Translational.Components.Mass mass(v(fixed = true, start = 0), m = 1300) annotation(
        Placement(transformation(extent = {{54, 42}, {74, 62}})));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, S = 2.2, Cx = 0.26, m = mass.m) annotation(
        Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 90, origin = {89, 29})));
      wbEHPTlib.MapBased.Partial.IceConnP ice annotation(
        Placement(transformation(extent = {{-98, 46}, {-78, 66}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1Conn bat(ECellMin = 0.9, ECellMax = 1.45, R0Cell = 0.0003, ns = 168, QCellNom = 2 * 6.5 * 3600.0, SOCInit = 0.6, ICellMax = 1e5, iCellEfficiency = 15 * 6.5) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-16, 0})));
      wbEHPTlib.SupportModels.ConnectorRelated.Conn d annotation(
        Placement(visible = true, transformation(extent = {{2, -40}, {28, -16}}, rotation = 0), iconTransformation(extent = {{4, -52}, {30, -28}}, rotation = 0)));
      wbEHPTlib.MapBased.ECUs.Ecu1 ECU(genLoopGain = 1.0) annotation(
        Placement(visible = true, transformation(origin = {-10, -41}, extent = {{-10, -9}, {10, 9}}, rotation = 0)));
      wbEHPTlib.MapBased.TwoFlangeConn mot( effTableName = "motEffTable", mapsFileName = "PSDmaps.txt",mapsOnFile = true) annotation(
        Placement(visible = true, transformation(extent = {{-28, 62}, {-8, 42}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(origin = {38, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {10, 26}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
      wbEHPTlib.MapBased.OneFlangeConn gen annotation(
        Placement(visible = true, transformation(extent = {{-38, 10}, {-58, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", k = 1, yMax = 1.8) annotation(
        Placement(visible = true, transformation(extent = {{-60, -50}, {-40, -30}}, rotation = 0)));
    equation
      connect(ECU.tauRef, driver.tauRef) annotation(
        Line(points = {{-22, -41}, {-29, -41}, {-29, -40}, {-39, -40}}, color = {0, 0, 127}));
      connect(carVel.v, driver.V) annotation(
        Line(points = {{78, -23}, {78, -58}, {-50, -58}, {-50, -51.2}}, color = {0, 0, 127}));
      connect(mot.conn1, ECU.conn1) annotation(
        Line(points = {{-27.2, 59.8}, {-27.2, 72}, {50, 72}, {50, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(gen.pin_p, bat.p) annotation(
        Line(points = {{-38, 24}, {-24, 24}, {-24, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(gen.pin_n, bat.n) annotation(
        Line(points = {{-38, 16}, {-8.5, 16}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(gen.flange_a, PSD.sun) annotation(
        Line(points = {{-58, 20}, {-70, 20}, {-70, 52}, {-60, 52}}));
      connect(gen.conn, ECU.conn1) annotation(
        Line(points = {{-58, 12.2}, {-58, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ground.p, bat.n) annotation(
        Line(points = {{0, 26}, {-8.5, 26}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(wheel.flangeT, mass.flange_a) annotation(
        Line(points = {{48, 52}, {54, 52}}, color = {0, 127, 0}));
      connect(wheel.flangeR, idealGear.flange_b) annotation(
        Line(points = {{28, 52}, {22, 52}}));
      connect(PSD.ring, mot.flange_a) annotation(
        Line(points = {{-40, 52}, {-34, 52}, {-28, 52}}));
      connect(idealGear.flange_a, mot.flange_b) annotation(
        Line(points = {{2, 52}, {-4, 52}, {-4, 52.2}, {-8, 52.2}}));
      connect(mot.pin_p, bat.p) annotation(
        Line(points = {{-22, 42.2}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(mot.pin_n, bat.n) annotation(
        Line(points = {{-14, 42}, {-14, 10.1}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(bat.conn, ECU.conn1) annotation(
        Line(points = {{-15.75, -10}, {-16, -10}, {-16, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ice.conn, ECU.conn1) annotation(
        Line(points = {{-88, 45.8}, {-88, 45.8}, {-88, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ECU.conn1, d) annotation(
        Line(points = {{-10, -32.18}, {-10, -28}, {15, -28}}, color = {255, 204, 51}, thickness = 0.5));
      der(EbatDel) = (bat.p.v - bat.n.v) * bat.n.i;
      der(EgenDelM) = gen.pin_p.i * (gen.pin_p.v - gen.pin_n.v) + gen.flange_a.tau * der(gen.flange_a.phi);
      der(Eroad) = dragForce.flange.f * der(dragForce.flange.s);
      der(EiceDel) = -ice.flange_a.tau * der(ice.flange_a.phi);
      der(Emot) = mot.flange_a.tau * der(mot.flange_a.phi) + mot.flange_b.tau * der(mot.flange_b.phi);
      Emass = 0.5 * mass.m * der(mass.flange_a.s) ^ 2;
      connect(PSD.carrier, ice.flange_a) annotation(
        Line(points = {{-60, 56}, {-70, 56}, {-70, 56}, {-78, 56}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(dragForce.flange, mass.flange_b) annotation(
        Line(points = {{89, 38}, {90, 38}, {90, 52}, {74, 52}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(carVel.flange, mass.flange_b) annotation(
        Line(points = {{78, -2}, {78, 52}, {74, 52}}, color = {0, 127, 0}, smooth = Smooth.None));
      annotation(
        __Dymola_experimentSetupOutput,
        Documentation(info = "<html>
<p>This model simulates a PSD based power train with a simple control logic, in its ECU:</p>
<p>It tries to make the ICE deliver the average load power, at its optimal speed.</p>
<p>This has two main inconveniences:</p>
<ul>
<li>the battery SOC is not controlled and tends to drift</li>
<li>in urban environments the power is too low to allow efficient drive without shutting off the engine.</li>
</ul>
<p>These inconveniences are addressed in subsequent models PSEcu2 and PSecu3 (see their infos).</p>
</html>"),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 2.8));
    end PSecu1;

    model PSecu2 "Full Power Split Device power train using Map-Based components"
      import Modelica.Constants.*;
      extends Modelica.Icons.Example;
      parameter Modelica.SIunits.AngularVelocity wIceStart = 50;
      Modelica.SIunits.Energy EbatDel "energy delivered by the battery";
      Modelica.SIunits.Energy EgenDelM "energy delivered by gen trough mechanical flange";
      Modelica.SIunits.Energy Eroad "mechanical energy absorbed by roas (friction & air)";
      Modelica.SIunits.Energy EiceDel "mechanical energy delivered by ice";
      Modelica.SIunits.Energy Emot, Emass;
      Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-50, 52})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
        Placement(transformation(extent = {{2, 42}, {22, 62}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor carVel annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {78, -12})));
      Modelica.Mechanics.Translational.Components.Mass mass(v(fixed = true, start = 0), m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{54, 42}, {74, 62}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, S = 2.2, Cx = 0.26, m = mass.m) annotation(
        Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 90, origin = {89, 29})));
      wbEHPTlib.MapBased.Partial.IceConnP ice(wIceStart = wIceStart, mapsFileName = "PSDmaps.txt") annotation(
        Placement(transformation(extent = {{-98, 46}, {-78, 66}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1Conn bat(ECellMin = 0.9, ECellMax = 1.45, R0Cell = 0.0003, ns = 168, QCellNom = 2 * 6.5 * 3600.0, SOCInit = 0.6, ICellMax = 1e5, iCellEfficiency = 15 * 6.5) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-16, 0})));
      wbEHPTlib.SupportModels.ConnectorRelated.Conn d annotation(
        Placement(visible = true, transformation(extent = {{2, -40}, {28, -16}}, rotation = 0), iconTransformation(extent = {{4, -52}, {30, -28}}, rotation = 0)));
      wbEHPTlib.MapBased.ECUs.Ecu2 ECU(socLoopGain = 1e5, genLoopGain = 1.0) annotation(
        Placement(visible = true, transformation(origin = {-10, -41}, extent = {{-10, -9}, {10, 9}}, rotation = 0)));
      wbEHPTlib.MapBased.TwoFlangeConn mot annotation(
        Placement(visible = true, transformation(extent = {{-28, 62}, {-8, 42}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(origin = {38, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {10, 26}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
      wbEHPTlib.MapBased.OneFlangeConn gen(mapsOnFile = true, mapsFileName = "PSDmaps.txt", effTableName = "genEffTable") annotation(
        Placement(visible = true, transformation(extent = {{-38, 10}, {-58, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver( CycleFileName = "NEDC.txt", k = 1,yMax = 1.8) annotation(
        Placement(visible = true, transformation(extent = {{-60, -50}, {-40, -30}}, rotation = 0)));
    equation
      connect(carVel.flange, mass.flange_b) annotation(
        Line(points = {{78, -2}, {78, 52}, {74, 52}}, color = {0, 127, 0}));
      connect(dragForce.flange, mass.flange_b) annotation(
        Line(points = {{89, 38}, {90, 38}, {90, 52}, {74, 52}}, color = {0, 127, 0}));
      connect(wheel.flangeT, mass.flange_a) annotation(
        Line(points = {{48, 52}, {54, 52}}, color = {0, 127, 0}));
      connect(ECU.tauReference, driver.tauRef) annotation(
        Line(points = {{-22, -41}, {-29, -41}, {-29, -40}, {-39, -40}}, color = {0, 0, 127}));
      connect(carVel.v, driver.V) annotation(
        Line(points = {{78, -23}, {78, -58}, {-50, -58}, {-50, -51.2}}, color = {0, 0, 127}));
      connect(mot.conn1, ECU.conn1) annotation(
        Line(points = {{-27.2, 59.8}, {-27.2, 72}, {50, 72}, {50, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(gen.pin_n, bat.n) annotation(
        Line(points = {{-38, 16}, {-8.5, 16}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(gen.flange_a, PSD.sun) annotation(
        Line(points = {{-58, 20}, {-70, 20}, {-70, 52}, {-60, 52}}));
      connect(gen.conn, ECU.conn1) annotation(
        Line(points = {{-58, 12.2}, {-58, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ground.p, bat.n) annotation(
        Line(points = {{0, 26}, {-8.5, 26}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(wheel.flangeR, idealGear.flange_b) annotation(
        Line(points = {{28, 52}, {22, 52}}));
      connect(PSD.ring, mot.flange_a) annotation(
        Line(points = {{-40, 52}, {-34, 52}, {-28, 52}}));
      connect(idealGear.flange_a, mot.flange_b) annotation(
        Line(points = {{2, 52}, {-4, 52}, {-4, 52.2}, {-8, 52.2}}));
      connect(mot.pin_p, bat.p) annotation(
        Line(points = {{-22, 42.2}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(mot.pin_n, bat.n) annotation(
        Line(points = {{-14, 42}, {-14, 10.1}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(bat.conn, ECU.conn1) annotation(
        Line(points = {{-15.75, -10}, {-16, -10}, {-16, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ice.conn, ECU.conn1) annotation(
        Line(points = {{-88, 45.8}, {-88, 45.8}, {-88, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ECU.conn1, d) annotation(
        Line(points = {{-10, -32.18}, {-10, -28}, {15, -28}}, color = {255, 204, 51}, thickness = 0.5));
      der(EbatDel) = (bat.p.v - bat.n.v) * bat.n.i;
      der(EgenDelM) = gen.pin_p.i * (gen.pin_p.v - gen.pin_n.v) + gen.flange_a.tau * der(gen.flange_a.phi);
      der(Eroad) = dragForce.flange.f * der(dragForce.flange.s);
      der(EiceDel) = -ice.flange_a.tau * der(ice.flange_a.phi);
      der(Emot) = mot.flange_a.tau * der(mot.flange_a.phi) + mot.flange_b.tau * der(mot.flange_b.phi);
      Emass = 0.5 * mass.m * der(mass.flange_a.s) ^ 2;
      connect(PSD.carrier, ice.flange_a) annotation(
        Line(points = {{-60, 56}, {-78, 56}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(gen.pin_p, bat.p) annotation(
        Line(points = {{-38, 24}, {-38, 34}, {-22, 34}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      annotation(
        __Dymola_experimentSetupOutput,
        Documentation(info = "<html>
<p>This model simulates a PSD based power train with a simple control logic, in its ECU:</p>
<p>it tries to make the ICE deliver the average load power, at its optimal speed and with a loop in the ECU that compensates SOC drift (improvement of Ecu2 over Ecu1, thus of PSecu2 over PSecu1).</p>
</html>"),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 2.8));
    end PSecu2;

    model PSecu3 "Full Power Split Device power train using Map-Based components"
      import Modelica.Constants.*;
      extends Modelica.Icons.Example;
      parameter Modelica.SIunits.AngularVelocity wIceStart = 50;
      parameter Real factorDebug = 100;
      // rad/s
      Modelica.SIunits.Energy EbatDel "energy delivered by the battery";
      Modelica.SIunits.Energy EgenDelM "energy delivered by gen trough mechanical flange";
      Modelica.SIunits.Energy Eroad "mechanical energy absorbed by roas (friction & air)";
      Modelica.SIunits.Energy EiceDel "mechanical energy delivered by ice";
      Modelica.SIunits.Energy Emot, Emass;
      Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-50, 52})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
        Placement(transformation(extent = {{2, 42}, {22, 62}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor carVel annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {78, -12})));
      Modelica.Mechanics.Translational.Components.Mass mass(m = 1300, v(fixed = true, start = 0)) annotation(
        Placement(transformation(extent = {{54, 42}, {74, 62}})));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(Cx = 0.26, S = 2.2, fc = 0.014, m = mass.m, rho = 1.226) annotation(
        Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 90, origin = {89, 29})));
      wbEHPTlib.MapBased.IceConnPOO ice(mapsFileName = "PSDmaps.txt", wIceStart = wIceStart) annotation(
        Placement(transformation(extent = {{-98, 46}, {-78, 66}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1Conn bat(ECellMin = 0.9, ECellMax = 1.45, R0Cell = 0.0003, ns = 168, QCellNom = 2 * 6.5 * 3600.0, SOCInit = 0.6, ICellMax = 1e5, iCellEfficiency = 15 * 6.5) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-16, 0})));
      wbEHPTlib.SupportModels.ConnectorRelated.Conn d annotation(
        Placement(visible = true, transformation(extent = {{2, -40}, {28, -16}}, rotation = 0), iconTransformation(extent = {{4, -52}, {30, -28}}, rotation = 0)));
      wbEHPTlib.MapBased.ECUs.Ecu3 ECU(genLoopGain = 1.0, socRef = 0.65) annotation(
        Placement(visible = true, transformation(origin = {-10, -41}, extent = {{-10, -9}, {10, 9}}, rotation = 0)));
      wbEHPTlib.MapBased.TwoFlangeConn mot annotation(
        Placement(visible = true, transformation(extent = {{-28, 62}, {-8, 42}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(origin = {38, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      wbEHPTlib.MapBased.OneFlangeConn gen annotation(
        Placement(visible = true, transformation(extent = {{-38, 10}, {-58, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(CycleFileName = "NEDC.txt", k = 1, yMax = 1.8) annotation(
        Placement(visible = true, transformation(extent = {{-60, -50}, {-40, -30}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {12, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ground.p, bat.n) annotation(
        Line(points = {{12, 20}, {-8, 20}, {-8, 10}, {-8, 10.1}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(ECU.tauReference, driver.tauRef) annotation(
        Line(points = {{-22, -41}, {-29, -41}, {-29, -40}, {-39, -40}}, color = {0, 0, 127}));
      connect(carVel.v, driver.V) annotation(
        Line(points = {{78, -23}, {78, -58}, {-50, -58}, {-50, -51.2}}, color = {0, 0, 127}));
      connect(mot.conn1, ECU.conn) annotation(
        Line(points = {{-27.2, 59.8}, {-27.2, 72}, {50, 72}, {50, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(gen.pin_p, bat.p) annotation(
        Line(points = {{-38, 24}, {-22, 24}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(gen.pin_n, bat.n) annotation(
        Line(points = {{-38, 16}, {-8.5, 16}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(gen.flange_a, PSD.sun) annotation(
        Line(points = {{-58, 20}, {-70, 20}, {-70, 52}, {-60, 52}}));
      connect(gen.conn, ECU.conn) annotation(
        Line(points = {{-58, 12.2}, {-58, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(wheel.flangeT, mass.flange_a) annotation(
        Line(points = {{48, 52}, {54, 52}}, color = {0, 127, 0}));
      connect(wheel.flangeR, idealGear.flange_b) annotation(
        Line(points = {{28, 52}, {22, 52}}));
      connect(PSD.ring, mot.flange_a) annotation(
        Line(points = {{-40, 52}, {-34, 52}, {-28, 52}}));
      connect(idealGear.flange_a, mot.flange_b) annotation(
        Line(points = {{2, 52}, {-4, 52}, {-4, 52.2}, {-8, 52.2}}));
      connect(mot.pin_p, bat.p) annotation(
        Line(points = {{-22, 42.2}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(mot.pin_n, bat.n) annotation(
        Line(points = {{-14, 42}, {-14, 10.1}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(bat.conn, ECU.conn) annotation(
        Line(points = {{-15.75, -10}, {-16, -10}, {-16, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ice.conn, ECU.conn) annotation(
        Line(points = {{-88, 46.2}, {-88, 46.2}, {-88, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ECU.conn, d) annotation(
        Line(points = {{-10, -32.18}, {-10, -28}, {15, -28}}, color = {255, 204, 51}, thickness = 0.5));
      der(EbatDel) = (bat.p.v - bat.n.v) * bat.n.i;
      der(EgenDelM) = gen.pin_p.i * (gen.pin_p.v - gen.pin_n.v) + gen.flange_a.tau * der(gen.flange_a.phi);
      der(Eroad) = dragForce.flange.f * der(dragForce.flange.s);
      der(EiceDel) = -ice.flange_a.tau * der(ice.flange_a.phi);
      der(Emot) = mot.flange_a.tau * der(mot.flange_a.phi) + mot.flange_b.tau * der(mot.flange_b.phi);
      Emass = 0.5 * mass.m * der(mass.flange_a.s) ^ 2;
      connect(PSD.carrier, ice.flange_a) annotation(
        Line(points = {{-60, 56}, {-78, 56}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(dragForce.flange, mass.flange_b) annotation(
        Line(points = {{89, 38}, {90, 38}, {90, 52}, {74, 52}}, color = {0, 127, 0}, smooth = Smooth.None));
      connect(carVel.flange, mass.flange_b) annotation(
        Line(points = {{78, -2}, {78, 52}, {74, 52}}, color = {0, 127, 0}, smooth = Smooth.None));
      annotation(
        __Dymola_experimentSetupOutput,
        Documentation(info = "<html>
<p>This model simulates a PSD based power train with a simple control logic, in its ECU:</p>
<p>it tries to make the ICE deliver the average load power, at its optimal speed and:</p>
<ol>
<li>has a loop in the ECU that compensates SOC drift (improvement Ecu1, thus over PSecu1)</li>
<li>has an hysteresis control in the ECU that implements ON/OFF strategy to avoid too low powers (and correspondingly high consumption) in urban environments (improvement of Ecu3 over Ecu1 and Ecu2, thus of PSecu3 over PSecu1 and PSecu2)</li>
</ol>
</html>"),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 2.8));
    end PSecu3;

    model PSecu2PA "Full Power Split Device power train using Map-Based components"
      import Modelica.Constants.*;
      extends Modelica.Icons.Example;
      parameter Modelica.SIunits.AngularVelocity wIceStart = 50;
      Modelica.SIunits.Energy EbatDel "energy delivered by the battery";
      Modelica.SIunits.Energy EgenDelM "energy delivered by gen trough mechanical flange";
      Modelica.SIunits.Energy Eroad "mechanical energy absorbed by roas (friction & air)";
      Modelica.SIunits.Energy EiceDel "mechanical energy delivered by ice";
      Modelica.SIunits.Energy Emot, Emass;
      Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-50, 52})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
        Placement(transformation(extent = {{2, 42}, {22, 62}})));
      Modelica.Mechanics.Translational.Sensors.SpeedSensor carVel annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 270, origin = {78, -12})));
      Modelica.Mechanics.Translational.Components.Mass mass(v(fixed = true, start = 0), m = 1300) annotation(
        Placement(visible = true, transformation(extent = {{54, 42}, {74, 62}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, S = 2.2, Cx = 0.26, m = mass.m) annotation(
        Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 90, origin = {89, 29})));
      wbEHPTlib.MapBased.Partial.IceConnP ice(wIceStart = wIceStart, mapsFileName = "PSDmaps.txt") annotation(
        Placement(transformation(extent = {{-98, 46}, {-78, 66}})));
      wbEHPTlib.SupportModels.Miscellaneous.Batt1Conn bat(ECellMin = 0.9, ECellMax = 1.45, R0Cell = 0.0003, SOCInit = 0.6, ICellMax = 1e5, iCellEfficiency = 15 * 6.5, ns = 100, QCellNom = 2 * 3600.0) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-16, 0})));
      wbEHPTlib.SupportModels.ConnectorRelated.Conn d annotation(
        Placement(visible = true, transformation(extent = {{2, -40}, {28, -16}}, rotation = 0), iconTransformation(extent = {{4, -52}, {30, -28}}, rotation = 0)));
      wbEHPTlib.MapBased.ECUs.Ecu2 ECU(genLoopGain = 1.0, socLoopGain = 2e4, powFiltT = 10) annotation(
        Placement(visible = true, transformation(origin = {-10, -41}, extent = {{-10, -9}, {10, 9}}, rotation = 0)));
      wbEHPTlib.MapBased.TwoFlangeConn mot annotation(
        Placement(visible = true, transformation(extent = {{-28, 62}, {-8, 42}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
        Placement(visible = true, transformation(origin = {38, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {10, 26}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
      wbEHPTlib.MapBased.OneFlangeConn gen annotation(
        Placement(visible = true, transformation(extent = {{-38, 10}, {-58, 30}}, rotation = 0)));
      wbEHPTlib.SupportModels.Miscellaneous.PropDriver driver(yMax = 1.8, CycleFileName = "nedc.txt", k = 1) annotation(
        Placement(visible = true, transformation(extent = {{-60, -50}, {-40, -30}}, rotation = 0)));
    equation
      connect(carVel.flange, mass.flange_b) annotation(
        Line(points = {{78, -2}, {78, 52}, {74, 52}}, color = {0, 127, 0}));
      connect(dragForce.flange, mass.flange_b) annotation(
        Line(points = {{89, 38}, {90, 38}, {90, 52}, {74, 52}}, color = {0, 127, 0}));
      connect(wheel.flangeT, mass.flange_a) annotation(
        Line(points = {{48, 52}, {54, 52}}, color = {0, 127, 0}));
      connect(ECU.tauReference, driver.tauRef) annotation(
        Line(points = {{-22, -41}, {-29, -41}, {-29, -40}, {-39, -40}}, color = {0, 0, 127}));
      connect(carVel.v, driver.V) annotation(
        Line(points = {{78, -23}, {78, -58}, {-50, -58}, {-50, -51.2}}, color = {0, 0, 127}));
      connect(mot.conn1, ECU.conn1) annotation(
        Line(points = {{-27.2, 59.8}, {-27.2, 72}, {50, 72}, {50, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(gen.pin_n, bat.n) annotation(
        Line(points = {{-38, 16}, {-8.5, 16}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(gen.flange_a, PSD.sun) annotation(
        Line(points = {{-58, 20}, {-58, 20}, {-70, 20}, {-70, 52}, {-60, 52}}));
      connect(gen.conn, ECU.conn1) annotation(
        Line(points = {{-58, 12.2}, {-58, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ground.p, bat.n) annotation(
        Line(points = {{0, 26}, {-8.5, 26}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(wheel.flangeR, idealGear.flange_b) annotation(
        Line(points = {{28, 52}, {22, 52}}));
      connect(PSD.ring, mot.flange_a) annotation(
        Line(points = {{-40, 52}, {-34, 52}, {-28, 52}}));
      connect(idealGear.flange_a, mot.flange_b) annotation(
        Line(points = {{2, 52}, {-4, 52}, {-4, 52.2}, {-8, 52.2}}));
      connect(mot.pin_p, bat.p) annotation(
        Line(points = {{-22, 42.2}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      connect(mot.pin_n, bat.n) annotation(
        Line(points = {{-14, 42}, {-14, 10.1}, {-8.5, 10.1}}, color = {0, 0, 255}));
      connect(bat.conn, ECU.conn1) annotation(
        Line(points = {{-15.75, -10}, {-16, -10}, {-16, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ice.conn, ECU.conn1) annotation(
        Line(points = {{-88, 45.8}, {-88, 45.8}, {-88, -20}, {-10, -20}, {-10, -32.18}}, color = {255, 204, 51}, thickness = 0.5));
      connect(ECU.conn1, d) annotation(
        Line(points = {{-10, -32.18}, {-10, -28}, {15, -28}}, color = {255, 204, 51}, thickness = 0.5));
      der(EbatDel) = (bat.p.v - bat.n.v) * bat.n.i;
      der(EgenDelM) = gen.pin_p.i * (gen.pin_p.v - gen.pin_n.v) + gen.flange_a.tau * der(gen.flange_a.phi);
      der(Eroad) = dragForce.flange.f * der(dragForce.flange.s);
      der(EiceDel) = -ice.flange_a.tau * der(ice.flange_a.phi);
      der(Emot) = mot.flange_a.tau * der(mot.flange_a.phi) + mot.flange_b.tau * der(mot.flange_b.phi);
      Emass = 0.5 * mass.m * der(mass.flange_a.s) ^ 2;
      connect(PSD.carrier, ice.flange_a) annotation(
        Line(points = {{-60, 56}, {-78, 56}}, color = {0, 0, 0}, smooth = Smooth.None));
      connect(gen.pin_p, bat.p) annotation(
        Line(points = {{-38, 24}, {-38, 34}, {-22, 34}, {-22, 10}, {-23.75, 10}}, color = {0, 0, 255}));
      annotation(
        __Dymola_experimentSetupOutput,
        Documentation(info = "<html>
<p>Like PSecu2, but solves the Proposed Activity that requires a very small battery and near-CVT (continuosuly-variable Transmission) operation. </p>
</html>"),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
        experiment(StartTime = 0, StopTime = 1400, Tolerance = 0.0001, Interval = 2.8));
    end PSecu2PA;

    package BasicPT "Basic Power Trains"
      extends Modelica.Icons.ExamplesPackage;

      model PSBasic0 "Modello a due macchine elettriche ideali (pure inerzie)"
        import Modelica.Constants.*;
        parameter Real vMass = 1300;
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
          Placement(transformation(extent = {{-24, 26}, {-4, 46}})));
        Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
          Placement(transformation(extent = {{0, 28}, {16, 44}})));
        Modelica.Mechanics.Translational.Components.Mass mass(m = 1300, v(fixed = true, start = 33.333333333333, displayUnit = "km/h")) annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {28, 26})));
        wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, S = 2.2, Cx = 0.26, m = mass.m, rho(displayUnit = "kg/m3") = 1.226) annotation(
          Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 90, origin = {27, -23})));
        Modelica.Mechanics.Rotational.Components.Inertia mot(w(fixed = false, displayUnit = "rpm"), J = 0.59) annotation(
          Placement(transformation(extent = {{-50, 26}, {-30, 46}})));
        Modelica.Mechanics.Translational.Sensors.PowerSensor Presis annotation(
          Placement(visible = true, transformation(origin = {28, -2}, extent = {{-8, -8}, {8, 8}}, rotation = 270)));
      equation
        connect(dragForce.flange, Presis.flange_b) annotation(
          Line(points = {{27, -16}, {27, -10}, {28, -10}}, color = {0, 127, 0}));
        connect(Presis.flange_a, mass.flange_b) annotation(
          Line(points = {{28, 6}, {28, 16}}, color = {0, 127, 0}));
        connect(wheel.flangeR, idealGear.flange_b) annotation(
          Line(points = {{0, 36}, {-4, 36}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGear.flange_a, mot.flange_b) annotation(
          Line(points = {{-24, 36}, {-24, 38}, {-26, 38}, {-26, 36}, {-30, 36}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(wheel.flangeT, mass.flange_a) annotation(
          Line(points = {{16, 36}, {28, 36}}, color = {0, 127, 0}));
        annotation(
          __Dymola_experimentSetupOutput,
          Documentation(info = "<html>
<h4>newInst OK</h4>
<h4>Modello di base per mostrare le prime funzionalit&agrave; di un sistema basato su PSD.</h4>
<p><br><u>Obiettivo finale</u>: simulazione di principio di Toyota Prius </p>
<p><u>Obiettivo di questa simulazione</u>: </p>
<p>Analizzare i punti di lavoro del PSD in un transitorio in cuil il veicolo decelera a seguito solo delle resistenze al moto d 120 km/h a 0; l&apos;ICE &egrave; mantenuto a velocit&agrave; costante da uno specifico controllo in retroazione, l&apos;inerzia di gen &egrave; aumentata di un fattore 20. </p>
</html>", revisions = "<html><head></head><body>No NewInst</body></html>"),
          conversion(noneFromVersion = ""),
          Diagram(coordinateSystem(extent = {{-80, -40}, {60, 60}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          Icon(coordinateSystem(extent = {{-80, -40}, {60, 60}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          experiment(StartTime = 0, StopTime = 1000, Tolerance = 0.0001, Interval = 0.5));
      end PSBasic0;

      model PSBasic1 "Modello a due macchine elettriche ideali (pure inerzie)"
        import Modelica.Constants.*;
        parameter Real vMass = 1300;
        Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-8, 0})));
        Modelica.Mechanics.Rotational.Components.Inertia ICE(J = 0.73, w(fixed = true, start = 300, displayUnit = "rad/s")) annotation(
          Placement(transformation(extent = {{-64, 14}, {-44, 34}})));
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
          Placement(transformation(extent = {{58, 0}, {78, 20}})));
        Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
          Placement(transformation(extent = {{82, 2}, {98, 18}})));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor vhVel annotation(
          Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 180, origin = {83, -21})));
        Modelica.Mechanics.Translational.Components.Mass mass(m = 1300, v(fixed = true, start = 33.333333333333, displayUnit = "km/h")) annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {110, 0})));
        wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, S = 2.2, Cx = 0.26, m = mass.m, rho(displayUnit = "kg/m3") = 1.226) annotation(
          Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 90, origin = {109, -49})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor iceW annotation(
          Placement(transformation(extent = {{-9, 9}, {9, -9}}, rotation = 180, origin = {1, 67})));
        Modelica.Mechanics.Rotational.Sources.Torque iceTau annotation(
          Placement(transformation(extent = {{-88, 14}, {-70, 32}})));
        Modelica.Blocks.Sources.Constant iceSetW(k = 300) annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 180, origin = {-2, 40})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor sunPow annotation(
          Placement(transformation(extent = {{-46, -18}, {-26, 2}})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor carrPow annotation(
          Placement(transformation(extent = {{-38, 16}, {-24, 30}})));
        Modelica.Mechanics.Rotational.Components.Inertia gen(w(fixed = false, displayUnit = "rpm"), J = 5) annotation(
          Placement(transformation(extent = {{-76, -18}, {-56, 2}})));
        Modelica.Mechanics.Rotational.Components.Inertia mot(w(fixed = false, displayUnit = "rpm"), J = 0.59) annotation(
          Placement(transformation(extent = {{32, 0}, {52, 20}})));
        Modelica.Blocks.Math.Feedback feedback annotation(
          Placement(transformation(extent = {{-22, 64}, {-42, 44}})));
        Modelica.Mechanics.Translational.Sensors.PowerSensor Presis annotation(
          Placement(visible = true, transformation(origin = {110, -28}, extent = {{-8, -8}, {8, 8}}, rotation = 270)));
        Modelica.Blocks.Math.Gain gain(k = 10) annotation(
          Placement(transformation(extent = {{-54, 44}, {-74, 64}})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor ringPow annotation(
          Placement(transformation(extent = {{24, -8}, {8, 8}})));
      equation
        connect(dragForce.flange, Presis.flange_b) annotation(
          Line(points = {{109, -42}, {109, -36}, {110, -36}}, color = {0, 127, 0}));
        connect(vhVel.flange, Presis.flange_a) annotation(
          Line(points = {{92, -21}, {92, -20}, {110, -20}}, color = {0, 127, 0}));
        connect(Presis.flange_a, mass.flange_b) annotation(
          Line(points = {{110, -20}, {110, -10}}, color = {0, 127, 0}));
        connect(wheel.flangeR, idealGear.flange_b) annotation(
          Line(points = {{82, 10}, {78, 10}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(iceTau.flange, ICE.flange_a) annotation(
          Line(points = {{-70, 23}, {-70, 24}, {-64, 24}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(sunPow.flange_b, PSD.sun) annotation(
          Line(points = {{-26, -8}, {-18, -8}, {-18, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(carrPow.flange_a, ICE.flange_b) annotation(
          Line(points = {{-38, 23}, {-42, 23}, {-42, 24}, {-44, 24}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(PSD.carrier, carrPow.flange_b) annotation(
          Line(points = {{-18, 4}, {-18, 23}, {-24, 23}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(carrPow.flange_b, iceW.flange) annotation(
          Line(points = {{-24, 23}, {18, 23}, {18, 67}, {10, 67}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(sunPow.flange_a, gen.flange_b) annotation(
          Line(points = {{-46, -8}, {-56, -8}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(idealGear.flange_a, mot.flange_b) annotation(
          Line(points = {{58, 10}, {58, 12}, {56, 12}, {56, 10}, {52, 10}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(iceSetW.y, feedback.u1) annotation(
          Line(points = {{-13, 40}, {-18, 40}, {-18, 54}, {-24, 54}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(iceW.w, feedback.u2) annotation(
          Line(points = {{-8.9, 67}, {-32, 67}, {-32, 62}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(feedback.y, gain.u) annotation(
          Line(points = {{-41, 54}, {-46, 54}, {-52, 54}}, color = {0, 0, 127}));
        connect(gain.y, iceTau.tau) annotation(
          Line(points = {{-75, 54}, {-86, 54}, {-96, 54}, {-96, 23}, {-89.8, 23}}, color = {0, 0, 127}));
        connect(wheel.flangeT, mass.flange_a) annotation(
          Line(points = {{98, 10}, {110, 10}}, color = {0, 127, 0}));
        connect(mot.flange_a, ringPow.flange_a) annotation(
          Line(points = {{32, 10}, {28, 10}, {28, 0}, {24, 0}}, color = {0, 0, 0}));
        connect(PSD.ring, ringPow.flange_b) annotation(
          Line(points = {{2, 0}, {8, 0}}, color = {0, 0, 0}));
        annotation(
          __Dymola_experimentSetupOutput,
          Documentation(info = "<html>
<h4>newInst OK</h4>
<h4>Modello di base per mostrare le prime funzionalit&agrave; di un sistema basato su PSD.</h4>
<p><br><u>Obiettivo finale</u>: simulazione di principio di Toyota Prius </p>
<p><u>Obiettivo di questa simulazione</u>: </p>
<p>Analizzare i punti di lavoro del PSD in un transitorio in cuil il veicolo decelera a seguito solo delle resistenze al moto d 120 km/h a 0; l&apos;ICE &egrave; mantenuto a velocit&agrave; costante da uno specifico controllo in retroazione, l&apos;inerzia di gen &egrave; aumentata di un fattore 20. </p>
</html>", revisions = "<html><head></head><body>No NewInst</body></html>"),
          conversion(noneFromVersion = ""),
          Diagram(coordinateSystem(extent = {{-100, -60}, {120, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          Icon(coordinateSystem(extent = {{-100, -60}, {120, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          experiment(StartTime = 0, StopTime = 1000, Tolerance = 0.0001, Interval = 0.5));
      end PSBasic1;

      model PSBasic2 "Modello a due macchine elettriche ideali (pure inerzie)"
        import Modelica.Constants.*;
        parameter Real vMass = 1300;
        parameter Real wICE = 167 "rad/s";
        Modelica.SIunits.Power genPow0 = genTau.tau * der(genTau.flange.phi);
        // rad/s
        Modelica.Blocks.Nonlinear.Limiter limTgen(uMax = 30) annotation(
          Placement(visible = true, transformation(origin = {-86, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Blocks.Math.Feedback feedback1 annotation(
          Placement(visible = true, transformation(extent = {{-22, -46}, {-42, -66}}, rotation = 0)));
        Modelica.Blocks.Math.Feedback feedback annotation(
          Placement(visible = true, transformation(extent = {{-14, 62}, {-34, 42}}, rotation = 0)));
        Modelica.Blocks.Math.Gain gain(k = -10) annotation(
          Placement(visible = true, transformation(extent = {{-54, -66}, {-74, -46}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sources.Torque genTau annotation(
          Placement(visible = true, transformation(extent = {{-82, -18}, {-64, 0}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.Inertia mot(w(fixed = false, displayUnit = "rpm"), J = 0.59) annotation(
          Placement(visible = true, transformation(extent = {{28, -8}, {48, 12}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.Inertia gen(w(fixed = false, displayUnit = "rpm"), J = 0.25) annotation(
          Placement(visible = true, transformation(extent = {{-58, -18}, {-38, 2}}, rotation = 0)));
        Modelica.Mechanics.Translational.Sensors.PowerSensor propPow annotation(
          Placement(visible = true, transformation(origin = {38, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor icePow annotation(
          Placement(visible = true, transformation(extent = {{-30, 12}, {-10, 32}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor genPow annotation(
          Placement(visible = true, transformation(extent = {{-32, -18}, {-12, 2}}, rotation = 0)));
        Modelica.Blocks.Sources.Constant iceSetW(k = 300) annotation(
          Placement(visible = true, transformation(origin = {18, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        Modelica.Mechanics.Rotational.Sources.Torque iceTau annotation(
          Placement(visible = true, transformation(extent = {{-80, 12}, {-60, 32}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor iceW annotation(
          Placement(visible = true, transformation(origin = {10, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        Modelica.Mechanics.Translational.Sensors.PowerSensor vhPow annotation(
          Placement(visible = true, transformation(origin = {93, -35}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(v(fixed = true, start = 0, displayUnit = "km/h"), m = 1300) annotation(
          Placement(visible = true, transformation(origin = {65, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor vhVel annotation(
          Placement(visible = true, transformation(origin = {70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
          Placement(visible = true, transformation(extent = {{84, -6}, {100, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
          Placement(visible = true, transformation(extent = {{56, -8}, {76, 12}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.Inertia ICE(J = 0.73, w(fixed = true, start = iceSetW.k, displayUnit = "rad/s")) annotation(
          Placement(visible = true, transformation(extent = {{-54, 12}, {-34, 32}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = 78 / 30) annotation(
          Placement(visible = true, transformation(origin = {10, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Sources.Constant vRef(k = 130 / 3.6) annotation(
          Placement(visible = true, transformation(origin = {-3, -55}, extent = {{-9, -9}, {9, 9}}, rotation = 180)));
        wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, S = 2.2, Cx = 0.26, m = mass.m) annotation(
          Placement(visible = true, transformation(origin = {107, -62}, extent = {{-9, -9}, {9, 9}}, rotation = 90)));
        Modelica.Blocks.Math.Gain gain1(k = 10) annotation(
          Placement(transformation(extent = {{-42, 42}, {-62, 62}})));
      equation
        connect(dragForce.flange, vhPow.flange_b) annotation(
          Line(points = {{107, -53}, {108, -53}, {108, -54}, {110, -54}, {110, -34}, {102, -34}, {102, -35}}, color = {0, 127, 0}));
        connect(vRef.y, feedback1.u1) annotation(
          Line(points = {{-12.9, -55}, {-11.9, -55}, {-11.9, -56}, {-24, -56}}, color = {0, 0, 127}));
        connect(PSD.ring, mot.flange_a) annotation(
          Line(points = {{20, 2}, {28, 2}}));
        connect(PSD.carrier, icePow.flange_b) annotation(
          Line(points = {{-5.55112e-16, 6}, {-5.55112e-16, 22}, {-10, 22}}));
        connect(genPow.flange_b, PSD.sun) annotation(
          Line(points = {{-12, -8}, {-12, -8}, {0, -8}, {0, 2}}));
        connect(icePow.flange_a, ICE.flange_b) annotation(
          Line(points = {{-30, 22}, {-34, 22}}));
        connect(iceTau.flange, ICE.flange_a) annotation(
          Line(points = {{-60, 22}, {-54, 22}}));
        connect(idealGear.flange_a, mot.flange_b) annotation(
          Line(points = {{56, 2}, {48, 2}}));
        connect(wheel.flangeR, idealGear.flange_b) annotation(
          Line(points = {{84, 2}, {76, 2}}));
        connect(wheel.flangeT, propPow.flange_a) annotation(
          Line(points = {{100, 2}, {112, 2}, {112, -22}, {20, -22}, {20, -34}, {28, -34}}, color = {0, 127, 0}));
        connect(vhVel.v, feedback1.u2) annotation(
          Line(points = {{59, -60}, {11.5, -60}, {11.5, -36}, {-31.25, -36}, {-31.25, -48}, {-32, -48}}, color = {0, 0, 127}));
        connect(vhVel.flange, vhPow.flange_a) annotation(
          Line(points = {{80, -60}, {80, -35}, {84, -35}}, color = {0, 127, 0}));
        connect(propPow.flange_b, mass.flange_a) annotation(
          Line(points = {{48, -34}, {54, -34}, {54, -35}}, color = {0, 127, 0}));
        connect(vhPow.flange_a, mass.flange_b) annotation(
          Line(points = {{84, -35}, {76, -35}}, color = {0, 127, 0}));
        connect(feedback.u2, iceW.w) annotation(
          Line(points = {{-24, 60}, {-24, 70}, {-1, 70}}, color = {0, 0, 127}));
        connect(icePow.flange_b, iceW.flange) annotation(
          Line(points = {{-10, 22}, {40, 22}, {40, 70}, {20, 70}}));
        connect(feedback.u1, iceSetW.y) annotation(
          Line(points = {{-16, 52}, {-8, 52}, {-8, 40}, {7, 40}}, color = {0, 0, 127}));
        connect(genPow.flange_a, gen.flange_b) annotation(
          Line(points = {{-32, -8}, {-38, -8}}));
        connect(gen.flange_a, genTau.flange) annotation(
          Line(points = {{-58, -8}, {-62, -8}, {-62, -9}, {-64, -9}}));
        connect(limTgen.y, genTau.tau) annotation(
          Line(points = {{-86, -29}, {-90, -29}, {-90, -10}, {-86, -10}, {-86, -9}, {-83.8, -9}}, color = {0, 0, 127}));
        connect(limTgen.u, gain.y) annotation(
          Line(points = {{-86, -52}, {-86, -56}, {-75, -56}}, color = {0, 0, 127}));
        connect(feedback1.y, gain.u) annotation(
          Line(points = {{-41, -56}, {-52, -56}}, color = {0, 0, 127}));
        connect(feedback.y, gain1.u) annotation(
          Line(points = {{-33, 52}, {-40, 52}}, color = {0, 0, 127}));
        connect(gain1.y, iceTau.tau) annotation(
          Line(points = {{-63, 52}, {-74, 52}, {-92, 52}, {-92, 22}, {-82, 22}}, color = {0, 0, 127}));
        annotation(
          experiment(StopTime = 120, Interval = 0.5),
          __Dymola_experimentSetupOutput,
          Documentation(info = "<html>
<h4>newInst OK</h4>
<h4>Modello di base per mostrare le prime funzionalit&agrave; di un sistema basato su PSD.</h4>
<p><br><u>Obiettivo finale</u>: simulazione di principio di Toyota Prius </p>
<p><u>Obiettivo di questa simulazione</u>: </p>
<p>Analizzare i punti di lavoro del PSD in un transitorio in cui il veicolo partendo da 0 accelera fino alla velocit&agrave; di 100 km/h, agendo sul gen; l&apos;ICE &egrave; mantenuto a velocit&agrave; costante da uno specifico controllo in retroazione. </p>
</html>"),
          conversion(noneFromVersion = ""),
          Icon(coordinateSystem(extent = {{-100, -80}, {120, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          Diagram(coordinateSystem(extent = {{-100, -80}, {120, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {18.4202, 31}, lineColor = {28, 108, 200}, extent = {{31.5798, 45}, {91.5798, 25}}, textString = "Valori iniziali:
velocità massa
velocità ICE
          "), Text(origin = {19.4729, 6}, lineColor = {238, 46, 47}, extent = {{30.5271, 54}, {88.5271, 30}}, textString = "Ripetere con valore più realistico 
di limTgen: 60Nm
Poi con velocità max di 120 km/h")}));
      end PSBasic2;

      model EleBalanceSim "Bilanciamento elettrico con simulazione (v. Epicicloidale.docx)"
        import Modelica.Constants.*;
        parameter Real vMass = 1300;
        parameter Real rho = 78 / 30;
        parameter Real sigma = 1 / rho;
        //  Real Tm = sigma * gear.flange_a.tau * wGen.w / (wMot.w + sigma * wGen.w);
        Real Tm = gear.flange_a.tau * ((1 + sigma) * wIce.w - wMot.w) / ((1 + sigma) * wIce.w);
        Real negPmot = -Tm * wMot.w;
        Modelica.Mechanics.Rotational.Components.IdealPlanetary PSD(ratio = rho) annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {-6, 0})));
        Modelica.Mechanics.Rotational.Components.Inertia ice(J = 0.73, w(fixed = true, start = 300, displayUnit = "rad/s")) annotation(
          Placement(transformation(extent = {{44, 50}, {64, 70}})));
        Modelica.Mechanics.Rotational.Components.IdealGear gear(ratio = 3.905) annotation(
          Placement(transformation(extent = {{46, -10}, {66, 10}})));
        Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
          Placement(transformation(extent = {{70, -8}, {86, 8}})));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor carVel annotation(
          Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 180, origin = {105, -43})));
        Modelica.Mechanics.Translational.Components.Mass mass(m = vMass, v(fixed = true, displayUnit = "km/h", start = 1.0)) annotation(
          Placement(transformation(extent = {{-9, -9}, {9, 9}}, rotation = 0, origin = {121, 1})));
        Modelica.Mechanics.Translational.Sensors.PowerSensor pResis annotation(
          Placement(transformation(extent = {{-6, -6}, {6, 6}}, rotation = 270, origin = {132, -24})));
        wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, m = vMass, S = 2.2, Cx = 0.26) annotation(
          Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 90, origin = {131, -49})));
        Modelica.Mechanics.Rotational.Sources.Torque tauIce annotation(
          Placement(transformation(extent = {{16, 50}, {36, 70}})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor pGen annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {-20, -28})));
        Modelica.Mechanics.Rotational.Sensors.PowerSensor pIce annotation(
          Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = 180, origin = {12, 36})));
        Modelica.Mechanics.Translational.Sensors.PowerSensor pProp annotation(
          Placement(transformation(extent = {{-7, -7}, {7, 7}}, rotation = 0, origin = {99, 1})));
        Modelica.Mechanics.Rotational.Components.Inertia gen(w(fixed = false, displayUnit = "rpm"), J = 0.1) annotation(
          Placement(transformation(extent = {{-46, -58}, {-26, -38}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor wGen annotation(
          Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = 0, origin = {-48, 0})));
        Modelica.Blocks.Math.Feedback fbIce annotation(
          Placement(transformation(extent = {{-46, 50}, {-26, 70}})));
        wbEHPTlib.SupportModels.MapBasedRelated.InertiaTq mot(J = 0.1) annotation(
          Placement(transformation(extent = {{12, -10}, {32, 10}})));
        Modelica.Mechanics.Rotational.Sources.Torque tauGen annotation(
          Placement(transformation(extent = {{-70, -58}, {-50, -38}})));
        Modelica.Blocks.Math.Gain gain(k = -50.0) annotation(
          Placement(transformation(extent = {{18, -76}, {-2, -56}})));
        Modelica.Blocks.Math.Feedback fbGen annotation(
          Placement(transformation(extent = {{50, -56}, {30, -76}})));
        Modelica.Blocks.Sources.Ramp ramp(height = 30, duration = 30) annotation(
          Placement(transformation(extent = {{84, -76}, {64, -56}})));
        Modelica.Blocks.Sources.RealExpression tauM(y = Tm) annotation(
          Placement(transformation(extent = {{50, -36}, {30, -16}})));
        Modelica.Blocks.Sources.Constant const(k = 300) annotation(
          Placement(transformation(extent = {{-72, 50}, {-52, 70}})));
        Modelica.Blocks.Math.Gain gain1(k = 100) annotation(
          Placement(transformation(extent = {{-18, 50}, {2, 70}})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor wIce annotation(
          Placement(transformation(extent = {{10, -10}, {-10, 10}}, rotation = -90, origin = {-36, 26})));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor wMot annotation(
          Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0, origin = {60, 20})));
      equation
        connect(wheel.flangeR, gear.flange_b) annotation(
          Line(points = {{70, 0}, {66, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(pResis.flange_a, mass.flange_b) annotation(
          Line(points = {{132, -18}, {132, 1}, {130, 1}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(carVel.flange, pResis.flange_a) annotation(
          Line(points = {{114, -43}, {114, -18}, {132, -18}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(dragForce.flange, pResis.flange_b) annotation(
          Line(points = {{131, -42}, {132, -42}, {132, -30}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(tauIce.flange, ice.flange_a) annotation(
          Line(points = {{36, 60}, {44, 60}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(pGen.flange_b, PSD.sun) annotation(
          Line(points = {{-20, -18}, {-20, 0}, {-16, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(pIce.flange_a, ice.flange_b) annotation(
          Line(points = {{22, 36}, {78, 36}, {78, 60}, {64, 60}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(PSD.carrier, pIce.flange_b) annotation(
          Line(points = {{-16, 4}, {-16, 36}, {2, 36}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(pProp.flange_a, wheel.flangeT) annotation(
          Line(points = {{92, 1}, {90, 1}, {90, 0}, {86, 0}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(pProp.flange_b, mass.flange_a) annotation(
          Line(points = {{106, 1}, {112, 1}}, color = {0, 127, 0}, smooth = Smooth.None));
        connect(pGen.flange_a, gen.flange_b) annotation(
          Line(points = {{-20, -38}, {-20, -48}, {-26, -48}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(wGen.flange, PSD.sun) annotation(
          Line(points = {{-38, 0}, {-16, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(gear.flange_a, mot.flange_b) annotation(
          Line(points = {{46, 0}, {32, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(PSD.ring, mot.flange_a) annotation(
          Line(points = {{4, 0}, {12, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(gen.flange_a, tauGen.flange) annotation(
          Line(points = {{-46, -48}, {-50, -48}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(fbGen.y, gain.u) annotation(
          Line(points = {{31, -66}, {20, -66}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(gain.y, tauGen.tau) annotation(
          Line(points = {{-3, -66}, {-78, -66}, {-78, -48}, {-72, -48}, {-72, -48}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(fbGen.u2, carVel.v) annotation(
          Line(points = {{40, -58}, {40, -43}, {95.1, -43}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(fbGen.u1, ramp.y) annotation(
          Line(points = {{48, -66}, {63, -66}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(tauM.y, mot.tau) annotation(
          Line(points = {{29, -26}, {16.55, -26}, {16.55, -10}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(const.y, fbIce.u1) annotation(
          Line(points = {{-51, 60}, {-44, 60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(wIce.flange, PSD.carrier) annotation(
          Line(points = {{-36, 16}, {-22, 16}, {-22, 4}, {-16, 4}}, color = {0, 0, 0}, smooth = Smooth.None));
        connect(wIce.w, fbIce.u2) annotation(
          Line(points = {{-36, 37}, {-36, 52}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(tauIce.tau, gain1.y) annotation(
          Line(points = {{14, 60}, {3, 60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(gain1.u, fbIce.y) annotation(
          Line(points = {{-20, 60}, {-27, 60}}, color = {0, 0, 127}, smooth = Smooth.None));
        connect(wMot.flange, mot.flange_b) annotation(
          Line(points = {{50, 20}, {38, 20}, {38, 0}, {32, 0}}, color = {0, 0, 0}, smooth = Smooth.None));
        annotation(
          experiment(StopTime = 40, Interval = 0.1),
          __Dymola_experimentSetupOutput,
          Documentation(info = "<html><head></head><body><p></p><h4>newInst OK</h4><h4>Modello di base per mostrare le prime funzionalità di un sistema basato su PSD.</h4><p></p>
    <p><br><u>Obiettivo finale</u>: simulazione di principio di Toyota Prius</p>
    <p><u>Obiettivo di questa simulazione</u>: </p>
    <p>Mostrare il controllo in bilanciamento elettrico.</p>
    </body></html>"),
          conversion(noneFromVersion = ""),
          Diagram(coordinateSystem(extent = {{-80, -80}, {140, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics),
          Icon(coordinateSystem(extent = {{-80, -80}, {140, 80}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})));
      end EleBalanceSim;

      model FreeBraking "Modello a due macchine elettriche ideali (pure inerzie)"
        import Modelica.Constants.*;
        parameter Real vMass = 1300;
        Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio = 3.905) annotation(
          Placement(visible = true, transformation(extent = {{-30, 10}, {-10, 30}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.IdealRollingWheel wheel(radius = 0.31) annotation(
          Placement(visible = true, transformation(extent = {{-6, 12}, {10, 28}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Components.Inertia mot(w(fixed = false, displayUnit = "rpm"), J = 0.59) annotation(
          Placement(visible = true, transformation(extent = {{-56, 10}, {-36, 30}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor wGen annotation(
          Placement(visible = true, transformation(origin = {-60, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
        Modelica.Mechanics.Translational.Sensors.PowerSensor Pprop annotation(
          Placement(visible = true, transformation(origin = {24, 20}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
        Modelica.Mechanics.Translational.Components.Mass mass(m = vMass, v(fixed = true, start = 33.333333333333, displayUnit = "km/h")) annotation(
          Placement(visible = true, transformation(origin = {46, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Translational.Sensors.PowerSensor Presis annotation(
          Placement(visible = true, transformation(origin = {62.5, -9.5}, extent = {{-9.5, -9.5}, {9.5, 9.5}}, rotation = 270)));
        wbEHPTlib.SupportModels.Miscellaneous.DragForce dragForce(fc = 0.014, rho = 1.226, m = vMass, S = 2.2, Cx = 0.26) annotation(
          Placement(visible = true, transformation(origin = {61.5, -36.5}, extent = {{-9.5, -9.5}, {9.5, 9.5}}, rotation = 90)));
        Modelica.Mechanics.Translational.Sensors.SpeedSensor Vvh annotation(
          Placement(visible = true, transformation(origin = {29, -17}, extent = {{-9, -9}, {9, 9}}, rotation = 180)));
      equation
        connect(Vvh.flange, Presis.flange_a) annotation(
          Line(points = {{38, -17}, {38, 0}, {62.5, 0}}, color = {0, 127, 0}));
        connect(dragForce.flange, Presis.flange_b) annotation(
          Line(points = {{61.5, -27}, {61.5, -19.5}, {62.5, -19.5}, {62.5, -19}}, color = {0, 127, 0}));
        connect(Presis.flange_a, mass.flange_b) annotation(
          Line(points = {{62.5, 0}, {62, 0}, {62, 20}, {56, 20}, {56, 20}, {56, 20}}, color = {0, 127, 0}));
        connect(dragForce.flange, Presis.flange_b) annotation(
          Line(points = {{61.5, -27}, {62.5, -27}, {62.5, -19}}, color = {0, 127, 0}));
        connect(Pprop.flange_b, mass.flange_a) annotation(
          Line(points = {{32, 20}, {36, 20}}, color = {0, 127, 0}));
        connect(Pprop.flange_a, wheel.flangeT) annotation(
          Line(points = {{16, 20}, {10, 20}}, color = {0, 127, 0}));
        connect(wGen.flange, mot.flange_a) annotation(
          Line(points = {{-60, -4}, {-60, 9}, {-60, 9}, {-60, 20}, {-58, 20}, {-58, 20}, {-56, 20}}));
        connect(idealGear.flange_a, mot.flange_b) annotation(
          Line(points = {{-30, 20}, {-36, 20}}));
        connect(wheel.flangeR, idealGear.flange_b) annotation(
          Line(points = {{-6, 20}, {-10, 20}}));
        annotation(
          __Dymola_experimentSetupOutput,
          Documentation(info = "<html><head></head><body><p></p><h4>Modello di base per mostrare le prime funzionalità di un sistema basato su PSD.</h4><p></p>
    <p><br><u>Obiettivo finale</u>: simulazione di principio di Toyota Prius</p>
    <p><u>Obiettivo di questa simulazione</u>: </p>
    <p>funzionamento con frenata naturale (=a seguito delle sole resistenze al moto) del veicolo di riferimento.</p>
    </body></html>"),
          conversion(noneFromVersion = ""),
          Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          Icon(coordinateSystem(extent = {{-100, -60}, {100, 60}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2})),
          experiment(StartTime = 0, StopTime = 300, Tolerance = 0.0001, Interval = 0.15));
      end FreeBraking;
    end BasicPT;
  end PSD;

  package VehicleData
    record Bus
      parameter Real m = 16000;
      parameter Real ratio = 6;
      parameter Real radius = 0.5715;
      parameter Real J = 5;
      //Costante H=5s
      parameter Real Cx = 0.65;
      parameter Real S = 6;
      parameter Real fc = 0.013;
      parameter Real rho = 1.226;
      parameter Real kContr = 1000;
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {162, 29, 33}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-94, 22}, {96, -16}}, lineColor = {162, 29, 33}, textString = "Bus"), Text(extent = {{-100, -104}, {100, -140}}, lineColor = {0, 0, 255}, textString = "%name")}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end Bus;

    record Car
      parameter Real m = 1300;
      parameter Real ratio = 3.905;
      parameter Real radius = 0.31;
      parameter Real J = 1.5;
      parameter Real Cx = 0.26;
      parameter Real S = 2.2;
      parameter Real fc = 0.014;
      parameter Real rho = 1.226;
      parameter Real kContr = 100;
      /* Per il calcolo di J abbiamo usato Tavv=2s, Pn=50kW, Wbase= quella 
           corrispondente a 36 km/h, quindi 252 rad/s*/
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {162, 29, 33}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-94, 22}, {96, -16}}, lineColor = {162, 29, 33}, textString = "Car"), Text(extent = {{-100, -104}, {100, -140}}, lineColor = {0, 0, 255}, textString = "%name")}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end Car;
  end VehicleData;
  annotation(
    uses(Modelica(version = "3.2.2")));
end wbEHVpkg;
