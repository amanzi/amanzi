import datetime

rcParams = {'font.family':'Times',
            'axes.2D.tickson':True,
            'axes.2D.title.fontscale':2.0,
            'axes.2D.label.fontscale':2.0,
            'axes.2D.x.title':"x-coordinate [m]",
            'axes.2D.y.title':"z-coordinate [m]",
            'axes.3D.tickson':False,
            'axes.3D.title.fontscale':2.0,
            'axes.3D.label.fontscale':2.0,
            'axes.3D.x.title':None,
            'axes.3D.y.title':None,
            'axes.3D.z.title':None,
            'legend.fontheight':0.05,
            'legend.scale':(1.,1.8),
            'legend.position':(0.02,0.76),
            'legend.title.fontheight':0.04,
            'legend.title.position':(0.018,0.82),
            'legend.minmax':True,
            'figsize':(2048,1536),
            'pseudocolor.linewidth':3,
            'contour.linewidth':1,
            'contour.color':'k',
            'time.format':'%b %d',
            'time.zero': datetime.datetime(year=2005, month=1, day=1),
            'time.location': (0.018,0.9),
            'time.fontheight':0.08,
            'time.round': None,
            'time.window': 1,
            'var.renames' : {"saturation liquid": "liquid saturation",
                             "saturation ice": "ice saturation",
                             "saturation gas": "gas saturation"
                             },
            'var.units' : {"saturation": "-",
                           "depth": "m",
                           "temperature": "K",
                           "pressure": "Pa"
                           },
            'var.limits' : {"saturation": (0.,1.),
                            "depth": (0.,0.01),
                            "temperature": (-25,10)
                            }
            }


rcParams_poster = {'legend.fontheight':0.025,
                   'legend.title.fontheight':0.025,
                   'legend.scale':(1.1,1.),
                   'time.fontheight':0.04,
                   }
            
import visit as v

# fonts
_fonts = {"Arial": 0,
         "Courier": 1,
         "Times": 2
         }
def getDefaultFont():
    return _fonts[rcParams['font.family']]


def getAnnotationAttributes():
    annot = v.AnnotationAttributes()
    annot.userInfoFlag = 0
    annot.databaseInfoFlag = 0

    # 3D
    annot.axes3D.triadFlag = 0
    annot.axes3D.bboxFlag = 0

    # clobber the names
    if rcParams['axes.3D.x.title'] is not None:
        annot.axes3D.xAxis.title.userTitle = 1
        annot.axes3D.xAxis.title.title = rcParams['axes.3D.x.title']
    else:
        annot.axes3D.xAxis.title.visible = 0

    if rcParams['axes.3D.y.title'] is not None:
        annot.axes3D.yAxis.title.userTitle = 1
        annot.axes3D.yAxis.title.title = rcParams['axes.3D.y.title']
    else:
        annot.axes3D.yAxis.title.visible = 0

    if rcParams['axes.3D.z.title'] is not None:
        annot.axes3D.zAxis.title.userTitle = 1
        annot.axes3D.zAxis.title.title = rcParams['axes.3D.z.title']
    else:
        annot.axes3D.zAxis.title.visible = 0


    # move the axes to outside edges
    annot.axes3D.tickLocation = annot.axes3D.OutsideEdges
    if not rcParams['axes.3D.tickson']:
        annot.axes3D.visible = 0

    # 2D
    if rcParams['axes.2D.x.title'] is not None:
        annot.axes2D.xAxis.title.userTitle = 1
        annot.axes2D.xAxis.title.title = rcParams['axes.2D.x.title']
    else:
        annot.axes2D.xAxis.title.visible = 0

    if rcParams['axes.2D.y.title'] is not None:
        annot.axes2D.yAxis.title.userTitle = 1
        annot.axes2D.yAxis.title.title = rcParams['axes.2D.y.title']
    else:
        annot.axes2D.xAxis.title.visible = 0

    if not rcParams['axes.2D.tickson']:
        annot.axes2D.visible = 0

    # Fonts
    fnum = getDefaultFont()
    annot.axes2D.xAxis.title.font.font = fnum
    annot.axes2D.xAxis.title.font.scale = rcParams['axes.2D.title.fontscale']
    annot.axes2D.xAxis.label.font.font = fnum
    annot.axes2D.xAxis.label.font.scale = rcParams['axes.2D.label.fontscale']
    annot.axes2D.yAxis.title.font.font = fnum
    annot.axes2D.yAxis.title.font.scale = rcParams['axes.2D.title.fontscale']
    annot.axes2D.yAxis.label.font.font = fnum
    annot.axes2D.yAxis.label.font.scale = rcParams['axes.2D.label.fontscale']
    annot.axes3D.xAxis.title.font.font = fnum
    annot.axes3D.xAxis.title.font.scale = rcParams['axes.3D.title.fontscale']
    annot.axes3D.xAxis.label.font.font = fnum
    annot.axes3D.xAxis.label.font.scale = rcParams['axes.3D.label.fontscale']
    annot.axes3D.yAxis.title.font.font = fnum
    annot.axes3D.yAxis.title.font.scale = rcParams['axes.3D.title.fontscale']
    annot.axes3D.yAxis.label.font.font = fnum
    annot.axes3D.yAxis.label.font.scale = rcParams['axes.3D.label.fontscale']
    annot.axes3D.zAxis.title.font.font = fnum
    annot.axes3D.zAxis.title.font.scale = rcParams['axes.3D.title.fontscale']
    annot.axes3D.zAxis.label.font.font = fnum
    annot.axes3D.zAxis.label.font.scale = rcParams['axes.3D.label.fontscale']
    
    return annot


def renameScalar(oldname):
    newname = oldname.split(".")[0].replace("_", " ").replace("-", " ")

    units = None
    try:
        newname = rcParams['var.renames'][newname]
    except KeyError:
        pass

    units = None
    try:
        units = rcParams['var.units'][newname]
    except KeyError:
        try:
            units = rcParams['var.units'][newname.split(" ")[-1]]
        except KeyError:
            pass

    if units is not None:
        newname = newname + " [%s]"%units
    return newname

def getLimits(name):
    newname = name.split(".")[0].replace("_", " ")
    try:
        return rcParams['var.limits'][newname]
    except KeyError:
        pass

    try:
        return rcParams['var.limits'][newname.split(" ")[-1]]
    except KeyError:
        pass

    try:
        return rcParams['var.limits'][newname.split(" ")[0]]
    except KeyError:
        pass

    return None
