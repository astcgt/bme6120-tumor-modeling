/** Rebuild the Traditional Chinese teaching deck with the bundled Artifact Tool.
 * Node runtime and Python may be overridden through RUNTIME_NODE_MODULES and
 * RUNTIME_PYTHON. No project Python or conda dependencies are changed.
 */
import fs from 'node:fs/promises';
import {spawnSync} from 'node:child_process';
import {teachModel} from './teach_model.mjs';
import {tellResults} from './tell_results.mjs';
import path from 'node:path';
import os from 'node:os';
import {pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';

const workspace = path.resolve(import.meta.dirname, '..', '..');
const runtime = path.join(os.homedir(), '.cache/codex-runtimes/codex-primary-runtime/dependencies');
const modules = process.env.RUNTIME_NODE_MODULES ?? path.join(runtime, 'node/node_modules');
process.env.RUNTIME_NODE_MODULES = modules;
const req = createRequire(path.join(modules, '__runtime__.cjs'));
const {Presentation, PresentationFile, FileBlob} = await import(pathToFileURL(req.resolve('@oai/artifact-tool')).href);
const skill = path.join(os.homedir(), '.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations');
const {finalizePresentation,applyPresentationChartFont} = await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const build = process.env.PRESENTATION_BUILD_DIR ?? path.join(workspace, '.presentation-build/code-v9');
const out = path.join(workspace, 'presentation/slides');
await fs.mkdir(build,{recursive:true});
await fs.mkdir(out,{recursive:true});

const INK='#173042', BLUE='#236B8E', RED='#A74F3D', GRAY='#657580', BG='#FCFBF8';
const FONT='Arial Unicode MS';
const deck=Presentation.create({slideSize:{width:1280,height:720}});
const records=[];
const nativeTables=[];
const chartOwners=[];
const results=JSON.parse(await fs.readFile(path.join(workspace,'code/presentation_build/data/results.json'),'utf8'));

function text(s,str,x,y,w,h,size=28,color=INK,bold=false,align='left'){
  const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  sh.text=str;
  sh.text.style={typeface:FONT,fontSize:size,color,bold,alignment:align,verticalAlignment:'top',autoFit:'none',wrap:'square',insets:{left:0,right:0,top:0,bottom:0}};
  return sh;
}
function slide(title,subtitle,notes,source='simulation.py: Config, simulate, mutate_offspring'){
  const s=deck.slides.add();s.background.fill=BG;
  const n=records.length+1;
  text(s,title,64,44,1152,72,43,INK,true);
  if(subtitle)text(s,subtitle,66,125,1148,66,26,GRAY);
  text(s,String(n).padStart(2,'0'),1170,683,48,26,18,GRAY,false,'right');
  s.speakerNotes.textFrame.setText(notes+'\n\n來源：本專案 '+source+'。參數為教學假設，未用臨床資料校準。');
  records.push({number:n,title,subtitle,notes,source});
  return s;
}
function eq(s,str,y=235,size=38,x=70,w=1140,h=75){return text(s,str,x,y,w,h,size,BLUE);}
function line(s,str,y,x=70,w=1140,size=28){return text(s,str,x,y,w,74,size);}
function foot(s,str){text(s,str,68,619,1100,60,22,GRAY);}
function node(s,str,x,y,w=260,h=90,color=BLUE){
  const sh=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'#F0F5F7',line:{fill:color,width:1.5}});
  sh.text=str;sh.text.style={typeface:FONT,fontSize:27,color:INK,alignment:'center',verticalAlignment:'middle',autoFit:'none',insets:{left:9,right:9,top:5,bottom:5}};
  return sh;
}
function connect(s,a,b,from='right',to='left',kind='straight'){
  return s.shapes.connect(a,b,{kind,fromSide:from,toSide:to,line:{fill:GRAY,width:2},tail:{type:'triangle',length:'med',width:'med'}});
}
function table(s,values,x,y,w,h,widths=null,size=27){
  const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,values:values.map(r=>r.map(String)),...(widths?{columnWidths:widths}:{})});
  nativeTables.push(records.length);
  t.styleOptions={headerRow:false,bandedRows:false};
  t.borders.assign({fill:'#D9E1E5',width:.6,style:'solid'});
  for(let i=0;i<values.length;i++)for(let j=0;j<values[0].length;j++){
    const c=t.getCell(i,j);c.fill=i===0?'#E8F0F4':BG;
    c.text.style={typeface:FONT,fontSize:size,color:INK,bold:i===0,alignment:'center',verticalAlignment:'middle',autoFit:'none',insets:{left:6,right:6,top:5,bottom:5}};
  }
  return t;
}
const P=[[1,0,0,0,0],[.3,0,.672,.028,0],[.09,0,.4032,.0168,.49],[.027,0,.18144,.00756,.784],[0,0,0,0,1]];
function vmul(v,m){return m[0].map((_,j)=>v.reduce((sum,z,i)=>sum+z*m[i][j],0));}
const pi0=[0,1,0,0,0],pi1=vmul(pi0,P),pi2=vmul(pi1,P);
for(const row of P)if(Math.abs(row.reduce((a,b)=>a+b,0)-1)>1e-12)throw Error('P row sum');
if(Math.abs(pi2[0]-.361236)>1e-12)throw Error('P squared');

function scientificChart(s,kind,categories,series,ytitle,{xTitle='世代',min=0,max=undefined,labels=false,yFormat='#,##0',majorUnit=undefined}={}){
 chartOwners.push(records.length);
 const ch=s.charts.add(kind,{
  position:{left:65,top:200,width:875,height:380},categories:categories.map(String),series,
  hasLegend:kind!=='bar',legend:{position:'bottom',textStyle:{typeface:FONT,fontSize:19,fill:INK}},
  chartFill:BG,chartLine:{fill:'none',width:0},plotAreaFill:BG,plotAreaLine:{fill:'none',width:0},
  xAxis:{title:{text:xTitle,textStyle:{typeface:FONT,fontSize:22,fill:INK}},tickLabelPosition:'low',textStyle:{typeface:FONT,fontSize:19,fill:INK},line:{fill:GRAY,width:1}},
  yAxis:{title:{text:ytitle,textStyle:{typeface:FONT,fontSize:22,fill:INK}},min,...(max?{max}:{}),...(majorUnit?{majorUnit}:{}),numberFormatCode:yFormat,
    textStyle:{typeface:FONT,fontSize:19,fill:INK},majorGridlines:{fill:'#E0E5E8',width:.6},line:{fill:GRAY,width:1}},
  lineOptions:{smooth:false},scatterOptions:{style:'marker'},
  barOptions:{direction:'column',grouping:'clustered',gapWidth:80},
  dataLabels:{showValue:labels,position:'outEnd',textStyle:{typeface:FONT,fontSize:22,fill:INK}},
 });
 applyPresentationChartFont(ch,{fontFamily:FONT});
 return ch;
}
function curve(name,values,color,percentile=false){return {name,values,fill:color,
 line:{fill:color,width:percentile?1.3:3,style:percentile?'dashed':'solid'},marker:{symbol:'none'}};}
function side(s,stat,body){text(s,stat,970,235,240,150,38,BLUE,true);text(s,body,970,408,240,180,24);}


const formulaAssets=[];
// 程式節錄：預設取自 simulation.py；範圍第三項可指定其他檔案（例如 experiments.py）。
// Code excerpts: taken from simulation.py by default; a third range item can name another file (e.g. experiments.py).
const sources={};
for(const file of ['simulation.py','experiments.py'])sources[file]=(await fs.readFile(path.join(workspace,'code',file),'utf8')).split('\n');
function snip(ranges){
 const lines=ranges.flatMap(([a,b,file='simulation.py'])=>sources[file].slice(a-1,b)).map(l=>l.split('#')[0].trimEnd()).filter(l=>l.trim()&&!l.trim().startsWith('def ')&&!l.trim().startsWith('"""'));
 const indentation=Math.min(...lines.map(l=>l.match(/^ */)[0].length));
 return lines.map(l=>l.slice(indentation)).join('\n');
}
function code(s,str,x,y,w,h){
 // Display long source lines across visual continuation lines. No executable code is changed.
 const rows=[];
 for(let raw of str.split('\n')){
  while(raw.length>55){
   let i=raw.lastIndexOf(' ',55);
   if(i<8){
    const breaks=[...raw.slice(0,55).matchAll(/[.(,]/g)].map(m=>m.index+1).filter(j=>j>=8);
    i=breaks.length?breaks.at(-1):55;
   }
   rows.push(raw.slice(0,i));raw='    '+raw.slice(i).trimStart();
  }
  rows.push(raw);
 }
 const size=Math.min(19, h/(rows.length*1.25));
 if(size<15.5)throw Error(`Code too dense on slide ${records.length}: ${rows.length} lines`);
 const sh=text(s,rows.join('\n'),x,y,w,h,size,INK);
 sh.text.style={...sh.text.style,typeface:'Menlo',fontSize:size,wrap:'none'};
}
function math(s,latex,x=70,y=235,w=1140,h=80){
 const id=formulaAssets.length;
 formulaAssets.push({slide:s,latex,file:path.join(build,'equations',`formula-${id}.png`),x,y,w,h});
 records.at(-1).notes+='\nLaTeX source: '+latex;
 s.speakerNotes.textFrame.setText(records.at(-1).notes+'\n\n來源：本專案 '+records.at(-1).source+'。');
}
teachModel({slide,text,line,foot,node,connect,table,math,code,snip,P,pi2,scientificChart,curve,results});
tellResults({slide,text,line,foot,table,math,scientificChart,curve,side,results});
{
 const s=slide('程式與數學如何共同解釋模擬結果','資料表示決定狀態，抽樣規則決定轉移，統計決定解讀',
 'population對應基因型計數向量，therapy_generation保存持續治療控制。growth_probabilities和il11_environment定義狀態相依環境，二項抽樣生成存活与子代，逐位元mutation等價於Q。這些規則共同定義巨大P，simulate只抽其中一條路徑。多次試驗估計結果分布；tests檢查實作語意，資料校準另需實驗證據。讀完應能追蹤程式每個變數的來源、分布與作用時點。',
 'simulation.py, experiments.py, test_simulation.py');
 table(s,[['程式物件／步驟','數學意義','可解釋的結果'],['population + therapy_generation','完整狀態 Sₜ','相同負荷的不同後續行為'],['survivors / triples / mutation','Binomial 抽樣與 Q','生長、選擇與新基因型'],['simulate + 多次試驗','P 的路徑抽樣','復發比例與試驗間變異'],['summarize + 取樣 + tests','指標、比較與不變量','可追溯並驗證的分析']],70,215,1140,340,[460,290,390],25);
 math(s,String.raw`\boldsymbol{\pi}_{t+1}=\boldsymbol{\pi}_tP,\qquad\mathrm{E}[\mathbf{N}_{t+1}\mid S_t=s]=\mathbf{n}M(s)`,70,565,1140,65);
}

await fs.writeFile(path.join(build,'equations.json'),JSON.stringify(formulaAssets.map(({latex,file})=>({latex,file})),null,2));
const rendered=spawnSync(path.join(workspace,'.venv/bin/python'),[path.join(workspace,'code/presentation_build/render_equations.py'),path.join(build,'equations.json')],{encoding:'utf8'});
if(rendered.status!==0)throw Error(rendered.stderr||rendered.stdout);
console.log(rendered.stdout.trim());
for(const f of formulaAssets){
 const bytes=await fs.readFile(f.file),aspect=bytes.readUInt32BE(16)/bytes.readUInt32BE(20);
 const width=Math.min(f.w,f.h*aspect),height=width/aspect;
 f.slide.images.add({blob:new Uint8Array(bytes),contentType:'image/png',alt:'LaTeX: '+f.latex,fit:'contain',position:{left:f.x+(f.w-width)/2,top:f.y+(f.h-height)/2,width,height}});
}
const candidate=path.join(build,'candidate.pptx');
await (await PresentationFile.exportPptx(deck)).save(candidate);
await fs.writeFile(path.join(build,'slides.json'),JSON.stringify(records,null,2));
const finalPath=path.join(out,process.env.PRESENTATION_FILENAME??'tumor_markov_code_results_zh_tw_v9.pptx');
const result=await finalizePresentation({workspaceDir:workspace,candidatePath:candidate,finalPath,
 pythonExecutable:process.env.RUNTIME_PYTHON??path.join(runtime,'python/bin/python3'),
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',
   ...[...new Set(nativeTables)].flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:records.length,requiredNativeTableOwnerSlides:[...new Set(nativeTables)],
 requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,nativeChartTargetApplication:'powerpoint',
 fontPolicy:{basis:'design',families:[FONT,'Menlo']},verifyArtifactToolImport:true,
 receiptPath:path.join(build,`${path.basename(finalPath)}.${Date.now()}.validation.json`)});
console.log(JSON.stringify({finalPath,result}));

// Render the checked final file, rather than using draft previews as proof.
const checked=await PresentationFile.importPptx(await FileBlob.load(finalPath));
for(let i=0;i<records.length;i++){
 const s=checked.slides.items[i];
 const blob=await checked.export({slide:s,format:'png',scale:1.25});
 await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
 const layout=await s.export({format:'layout'});
 await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.layout.json`),await layout.text());
}
console.log(`Rendered all ${records.length} final slides.`);
