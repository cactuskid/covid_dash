
import json
import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
from dash_bio_utils import circos_parser as cp
from dash.dependencies import Input, Output, State
import dash_daq as daq
from dash_bio_utils import pdb_parser as parser, styles_parser as sparser
import dash_bio

import pandas as pd
from itertools import combinations

import numpy as np
import os
import wget
import glob
import time

import pandas as pd
import tempfile
from shutil import copy2
from textwrap import dedent as s


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

with open( './circoslayout.json'  , 'r') as circosin :
    circos_graph_data = json.loads( circosin.read() , parse_int=int )

clusterdf = pd.read_csv( './gisaid_hcov-2020_07_30.QC.NSoutlier.filter.deMaiomask.EPIID.HF.noambig.aln_codon_clustersclusterpositions_columns.csv')


models = '5x58,6cs0,6sc1,6nb6, 2jw8, 2xab, 4aud1, 1yo4, 2acf, 2wct, 3vc8, 2gt7, 3ee7, 2g9t, 3ee7, 2g9t, 6jyt, 1ysy, 6nur, 2g9t, 5c8u, 2g9t, 2xyq, 4mm3, 6cs2, 6acg, 6acj, 6ack, 2dd8, 2ghw, 6nb6, 6nb7'
models = models.split(',')
templatedir = './templates/'
clear = False
if not os.path.isdir(templatedir):
    os.mkdir(templatedir)

if clear == True:
    files = glob.glob(templatedir+'*.pdb')
    for f in files:
        os.remove(f)
dl_url = 'http://files.rcsb.org/download/'
dl_url_err = 'http://files.rcsb.org/download/'
structs = {}
already = glob.glob( './templates/*.pdb' )

for m in models:
    structfile = './templates/'+m.upper().strip()+'.pdb'
    if structfile not in already:
        print(m)
        time.sleep(1)
        try:
            wget.download(url = dl_url + m.strip() +'.pdb' , out =structfile)
            structs[m] = structfile
        except:
            try:
                wget.download(url = dl_url + m.strip() +'.pdb' , out =structfile)
                structs[m] = structfile
            except:
                print('err', m )
    else:
        structs[m.strip()] = structfile

#open all the complexes and pull sequences
#convert to json
out = 'outannot.txt'
annotation = pd.read_csv( out , header = None )
annotation.columns = [ 'qseqid' , 'sseqid' , 'qlen' ,  'slen' , 'qstart' , 'qend' ,  'qframe' , 'evalue' ]
annotation = annotation[ annotation['evalue'] < 10**-3 ]
annotation = annotation[annotation.evalue < 10 ** -20 ]
annotation = annotation.drop_duplicates(['sseqid'])
annotation['struct'] = annotation.sseqid.map(lambda x : x.split('|')[0])
annotation['chain'] = annotation.sseqid.map(lambda x : x.split('|')[1])

def files_data_style(content):
    fdat = tempfile.NamedTemporaryFile(suffix=".js", delete=False, mode='w+')
    fdat.write(content)
    dataFile = fdat.name
    fdat.close()
    return dataFile


# Create the model data from the decoded contents
fname = structs[list(structs.keys())[3]]
modata = parser.create_data(fname)
fmodel = files_data_style(modata)
with open(fmodel) as fm:
    mdata = json.load(fm)

datstyle = sparser.create_style(fname, 'cartoon', 'residue')
fstyle = files_data_style(datstyle)
with open(fstyle) as sf:
    data_style = json.load(sf)

maxnum = 0
for serial in data_style:
    if int(serial)>maxnum :
        maxnum = int(serial)

maxserial = 0
for atom_id in mdata['atoms']:
    if int(atom_id['serial'])>maxserial :
        maxserial = int(atom_id['serial'])

circosfinal = {}
for key in circos_graph_data:
    try:
        int(key)
        circosfinal[int(key)] = circos_graph_data[key]
    except:
        circosfinal[key] = circos_graph_data[key]
circos_graph_data = circosfinal

app = dash.Dash(__name__, external_stylesheets=external_stylesheets )

server = app.server

data_info = { os.path.abspath( structs[s])  : {
        'name': 'covid prot'+str(i),
        'description': 'test description',
        'link':'test link' } for i,s in enumerate(structs) }
stacktrack = {
                    'type': 'STACK',
                    'data': circos_graph_data['stack'],
                    'config': {
                        'innerRadius': .9 ,
                        'outerRadius': 1.2,
                        'thickness': 4,
                        'margin': 0,
                        'direction': 'out',
                        'strokeWidth': 0,
                        'opacity': 0.5,
                        'tooltipContent': {'name': 'chr'}
                        ,'color': {
                            'conditional': {
                                'end': 'end',
                                'start': 'start',
                                'value': [ n for n in range( 30000 , 1000)  ],
                                'color': [
                                    'red',
                                    'black',
                                    '#fff',
                                    '#999',
                                    '#BBB',
                                ],
                            }

                        }
                }
}

structtrack = {
                    'type': 'STACK',
                    'data': circos_graph_data['struct_stack'],
                    'config': {
                        'innerRadius': 0.7,
                        'outerRadius': 0.9,
                        'thickness': 4,
                        'margin': 0,
                        'direction': 'out',
                        'strokeWidth': 0,
                        'opacity': 0.5,
                        'tooltipContent': {'name': 'chr'}
                        }
                }
def header_colors():
    return {
        'bg_color': '#e7625f',
        'font_color': 'white'
    }

def description():
    return 'Molecule visualization in 3D - perfect for viewing ' \
           'biomolecules such as proteins, DNA and RNA. Includes ' \
           'stick, cartoon, and sphere representations.'

app.layout = html.Div(
        id='mol3d-body',
        className = 'row',
        children=[

######################### circos control pannel #####################################
        html.Div([

        html.Div([
        dashbio.Circos(
        id='my-dashbio-circos',
        layout= circos_graph_data['layout'],
        selectEvent={"0": "hover", "1": "click", "2": "both"},
        tracks=[

            {
            'type': 'CHORDS',
            'data': circos_graph_data[0],

            'config': {
                'tooltipContent': {
                    'source': 'source',
                    'sourceID': 'id',
                    'target': 'target',
                    'targetID': 'id',
                    'targetEnd': 'end'
                }
            }
        }
            ,
            stacktrack
            ,
            structtrack


        ] )

        ], id = 'circos' , className="six columns"  )



        ,

        dcc.Store(
                id='mol3d-color-storage',
                data={}
            ),

        html.Div(
        id='mol3d-biomolecule-viewer',
        children= [

############################# 3d viewer ####################################################
        dash_bio.Molecule3dViewer(
        id='mol-3d',
        selectionType='atom',
        modelData=mdata,
        styles = data_style ,
        selectedAtomIds=[],
        backgroundOpacity='0',
        atomLabelsShown=False,
        )


        ] , className="six columns")
        ,

####################### control pannel ############################################

        html.Div( className="six columns" , children = [

            html.Div(
                id='mol3d-control-tabs',
                className='control-tabs',
                children=[
                    dcc.Tabs(id='mol3d-tabs', value='what-is', children=[
                        dcc.Tab(
                            label='cluster select',
                            value='what-is',
                            children=html.Div(className='control-tab', children=[
                            "Graph type:",
                            dcc.Dropdown(
                                id='histogram-chords',
                                options=[
                                    {'label': x, 'value': x}
                                    for x in list(circos_graph_data.keys()) if type(x) is int
                                    ],
                                value=0
                                        ),

                            "Event data:",
                            html.Div(id='circos-output')

                            ] )
                        )

                        ,

                        dcc.Tab(
                            label='Structure select',
                            value='upload-download',
                            children=html.Div(className=' six columns', children=[
                                html.Div(
                                    title='Select molecule to view',
                                    className="app-controls-block",
                                    children=[
                                        html.Div(className='app-controls-name',
                                                 children='Select structure'),
                                        dcc.Dropdown(
                                            id='dropdown-demostr',
                                            options=[ { 'label':s,  'value':os.path.abspath(structs[s]) } for s in structs ],
                                            value= os.path.abspath(structs[list(structs.keys())[3]])
                                        ),
                                    ],
                                ),
                            html.Div(id='mol3d-data-info')
                            ])
                        ),

                        dcc.Tab(
                            label='Structure Info and coloring',
                            value='view-options',
                            children=html.Div(className='control-tab', children=[
                                # Textarea container to display the selected atoms
                                html.Div(
                                    title='view information about selected atoms '
                                    'of biomolecule',
                                    className="app-controls-block",
                                    id="mol3d-selection-display",
                                    children=[
                                        html.P(
                                            "Selection",
                                            style={
                                                'font-weight': 'bold',
                                                'margin-bottom': '10px'
                                            }
                                        ),
                                        html.Div(id='mol3d-selection-output'),
                                    ]
                                ),
                                # Dropdown to select chain representation
                                # (sticks, cartoon, sphere)
                                html.Div(
                                    title='select style for molecule representation',
                                    className="app-controls-block",
                                    id='mol3d-style',
                                    children=[
                                        html.P(
                                            'Style',
                                            style={
                                                'font-weight': 'bold',
                                                'margin-bottom': '10px'
                                            }
                                        ),
                                        dcc.Dropdown(
                                            id='dropdown-styles',
                                            options=[
                                                {'label': 'Sticks', 'value': 'stick'},
                                                {'label': 'Cartoon', 'value': 'cartoon'},
                                                {'label': 'Spheres', 'value': 'sphere'},
                                            ],
                                            value='cartoon'
                                        ),
                                    ],
                                ),

                                # Dropdown to select color of representation
                                html.Div(
                                    title='select color scheme for viewing biomolecule',
                                    className="app-controls-block",
                                    id='mol3d-style-color',
                                    children=[
                                        html.P(
                                            'Color',
                                            style={
                                                'font-weight': 'bold',
                                                'margin-bottom': '10px'
                                            }
                                        ),
                                        dcc.Dropdown(
                                            id='dropdown-style-color',
                                            options=[
                                                {'label': 'Atom',
                                                 'value': 'atom'},
                                                {'label': 'Residue identity',
                                                 'value': 'residue'},
                                                {'label': 'Residue type',
                                                 'value': 'residue_type'},
                                                {'label': 'Chain',
                                                 'value': 'chain'},
                                            ],
                                            value='residue'
                                        ),
                                        dcc.Dropdown(
                                            id='mol3d-coloring-key',
                                            options=[]
                                        ),

                                    ],
                                ),
                                html.Div(
                                    title='Customize molecule coloring.',
                                    className="app-controls-block",
                                    children=[
                                        html.P(
                                            id='mol3d-customize-coloring',
                                            style={
                                                'font-weight': 'bold',
                                                'margin-bottom': '10px'
                                            }
                                        ),
                                        daq.ColorPicker(
                                            id='mol3d-coloring-value',
                                            size=315
                                        ),
                                    ]
                                ),
                            ]),
                        ),

                    ]),
                ]  ),
        ])
        ],className="row") , dcc.Store(id='structswitch' )   ])



########################circos functions ########################
@app.callback(

    Output('circos-output', 'children'),
    [Input('my-dashbio-circos', 'eventDatum')]

)
def update_output(value):
    if value is not None:
        return [html.Div('{}: {}'.format(v.title(), value[v]))
                for v in value.keys()]
    return 'There are no event data. Click or hover on a data point to get more information.'

####change the circos content
@app.callback(
    Output('my-dashbio-circos', 'tracks'),
    [Input('histogram-chords', 'value')],
    state=[State('my-dashbio-circos', 'tracks')]
)
def change_graph_type(value, current):
    current[0].update(
            {

            'type': 'CHORDS',
            'data': circos_graph_data[int(value)] ,
            'config': {
                'tooltipContent': {
                    'source': 'source',
                    'sourceID': 'id',
                    'target': 'target',
                    'targetID': 'id',
                    'targetEnd': 'end'
                }
            }
        }
    )
    current[1].update(stacktrack)
    current[2].update(structtrack)
    return current


######### Function to create the modelData and style files for molecule visualization
@app.callback(
    Output('mol3d-data-info', 'children'),
    [Input('dropdown-demostr', 'value')]
)
def show_data(molecule_selected):
    if molecule_selected in data_info.keys():
        mol = data_info[molecule_selected]
        return [
            html.H4(mol['name']),
            mol['description'],
            html.A(
                '(source)',
                href=mol['link']
            )
        ]
    return ''

# Callback for updating dropdown options
@app.callback(
    Output('mol3d-coloring-key', 'options'),
    [Input('dropdown-style-color', 'value')]
)
def update_color_options(mol_style):
    color_dict_keys = {
        'atom': list(sparser.ATOM_COLOR_DICT.keys()),
        'residue': list(sparser.RESIDUE_COLOR_DICT.keys()),
        'residue_type': list(sparser.RESIDUE_TYPE_COLOR_DICT.keys()),
        'chain': list(sparser.CHAIN_COLOR_DICT.keys())
    }
    options = [{'label': k.upper(), 'value': k}
               for k in color_dict_keys[mol_style]]
    return options

@app.callback(
    [Output('structswitch','children'),Output("structswitch", "data")],
    [ Input('dropdown-demostr', 'value'),Input('histogram-chords', 'value') , Input('mol-3d' , 'modelData') ] ,
)
def set_color_cluster( structfile, cluster , mdata ):
    #blue and yellow high contrast
    #find match in annotation and blast ouput
    #geno positions from cluster file
    positions = [ int(p) for p in clusterdf.loc[cluster] if ~np.isnan(p) ]
    #select struct
    print(structfile)
    struct = structfile.split('/')[-1].replace('.pdb','').lower()
    geno_annot = annotation[ struct == annotation.struct ]
    print(geno_annot)
    residues = {}
    for i,r in geno_annot.iterrows():
        if r.chain not in residues:
            residues[r.chain] = []
        for p in positions:
            if p < r.qend and p > r.qstart:
                print(p)
                #match codon to amino acid position
                #todo: triple check this!!!!
                residues[r.chain].append( int( (p - r.qstart) / 3 ) )


    currentselection = []
    chaindiffs = {}
    chaincount = 0

    for atom_data in mdata['atoms']:
        for chain in residues:
            ###todo : triple check this
            if atom_data['chain'] == chain  and int(atom_data['residue_name'][3:]) in residues[chain]:
                currentselection.append(atom_data['serial'])
            #color select
    return structfile , currentselection

# Callback for molecule visualization based on uploaded PDB file
@app.callback(
    Output('mol3d-biomolecule-viewer', 'children'),
    [Input("structswitch", "data") , Input('dropdown-demostr', 'value'),
    Input('dropdown-styles', 'value'),
    Input('dropdown-style-color', 'value'),
    Input('mol3d-color-storage', 'modified_timestamp')],
    [State('mol3d-color-storage', 'data') ],
)
def use_upload( selection , demostr, mol_style, color_style, mt, custom_colors  ):
    print(selection)
    fname = demostr

    copy2(demostr, './str.pdb')
    fname = './str.pdb'
    # Create the model data from the decoded contents

    print(custom_colors)
    modata = parser.create_data(fname)
    fmodel = files_data_style(modata)
    with open(fmodel) as fm:
        mdata = json.load(fm)
    #all blue
    datstyle = sparser.create_style(fname, mol_style, color_style, **custom_colors)
    fstyle = files_data_style(datstyle)
    with open(fstyle) as sf:
        data_style = json.load(sf)

    print('selection finish' , selection)
    data_style = [ {'color': 'blue', 'visualization_type': 'cartoon'}  if i not in selection else  {'color': 'yellow', 'visualization_type': 'cartoon'}  for i,atom in enumerate(data_style) ]

    # Return the new molecule visualization container
    return dash_bio.Molecule3dViewer(
        id='mol-3d',
        selectionType='atom',
        modelData=mdata,
        styles=data_style,
        selectedAtomIds= [] ,
        backgroundOpacity='0',
        atomLabelsShown=False,
        )

#### add struct select callback from circos#######
# Callback to print details of each selected atom of the biomolecule
@app.callback(
    Output("mol3d-selection-output", "children"),
    [Input("mol-3d", "selectedAtomIds"),
     Input("mol-3d", "modelData")]
)
def selout(selected_atom_ids, model_data):
    residue_summary = []
    for atom_id in selected_atom_ids:
        res_info = model_data['atoms'][atom_id]
        residues = {
            "residue": res_info['residue_name'],
            "atom": res_info['name'],
            "chain": res_info['chain'],
            "xyz": res_info['positions']}

        residue_summary += [html.P('{}: {}'.format(
            key, str(residues[key]))) for key in residues]
        residue_summary.append(html.Br())
    if len(residue_summary) == 0:
        residue_summary.append("No atoms have been selected. Click \
        on an atom to select it.")
    return html.Div( residue_summary )

if __name__ == '__main__':
    app.run_server(debug = True )
