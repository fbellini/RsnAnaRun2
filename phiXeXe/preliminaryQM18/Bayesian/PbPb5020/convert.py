#!/usr/bin/python


from ROOT import TGraph, TCanvas
import os


def GetValues(fname, merge=[]):
    file = open(fname, "r")
    lines = file.readlines()
    lines = [w.replace('\n', '') for w in lines]
    values = []
    mcounter = 0
    for i in lines:
        j = i.split()
        print j
        for k in j[2:]:
            print k
            if len(merge) > 0 and mcounter < len(merge) and len(values) == merge[mcounter] + 1:
                mcounter += 1
                print "At bin {} summing {} to {}".format(len(values), k, values[-1])
                values[-1] += float(k)
                values[-1] /= 2.
            else:
                values.append(float(k))
    return values


def GetdNdy(part="pion"):
    print "dNdy for " + part
    return GetValues("dN_dy_{}.dat".format(part))


def GetMeanPt(part="pion"):
    print "MeanPt for " + part
    return GetValues("mean_pT_{}.dat".format(part))


def GetdNdeta():
    print "dNdeta"
    path = os.getcwd()
    if "PbPb2760" in path:
        return GetValues("dNch_deta.dat", [])
    else:
        return GetValues("dNch_deta.dat", [0, 1])


def FormGraph(x, y):
    print "x="
    print x
    print "y="
    print y
    g = TGraph(len(y))
    for i in enumerate(y):
        g.SetPoint(i[0], x[i[0]], i[1])
    g.Print("all")
    return g


def FormdNdy(part="pion"):
    dNdy = GetdNdy(part)
    dNdeta = GetdNdeta()
    return FormGraph(dNdeta, dNdy)


def FormMeanPt(part="pion"):
    MeanPt = GetMeanPt(part)
    dNdeta = GetdNdeta()
    return FormGraph(dNdeta, MeanPt)


def FormRatio(part=["proton", "pion"]):
    dNdy = [GetdNdy(part[0]), GetdNdy(part[1])]
    dNdeta = GetdNdeta()
    Ratio = []
    for i in enumerate(dNdy[0]):
        Ratio.append(i[1] / dNdy[1][i[0]])
    return FormGraph(dNdeta, Ratio)


def convert():
    L = "Pion Kaon Proton"
    L = L.split()
    for i in L:
        g = FormdNdy(i.lower())
        g.SaveAs("BayesianPrediction{}.root".format(i))
        g = FormMeanPt(i.lower())
        g.SaveAs("BayesianPredictionMeanPt{}.root".format(i))


def convertRatio():
    L = "Kaon Proton"
    L = L.split()
    for i in L:
        g = FormRatio([i.lower(), "pion"])
        g.SaveAs("BayesianPrediction{}_Over_Pion.root".format(i))


convert()
convertRatio()
