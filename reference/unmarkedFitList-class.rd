<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Class "unmarkedFitList" — unmarkedFitList-class • unmarked</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js" integrity="sha512-v2CJ7UaYy4JwqLDIrZUI/4hqeoQieOmAZNXBeQyjo21dadnwR+8ZaIJVT8EE2iyI61OV8e6M8PP2/4hpQINQ/g==" crossorigin="anonymous" referrerpolicy="no-referrer"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Class " unmarkedfitlist unmarkedfitlist-class><meta property="og:description" content="Class to hold multiple fitted models from one of
	unmarked's fitting functions"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">


    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">unmarked</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.4.1.9010</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li>
  <a href="../articles/unmarked.html">Get started</a>
</li>
<li>
  <a href="../reference/index.html">Reference</a>
</li>
<li class="dropdown">
  <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
    Articles

    <span class="caret"></span>
  </a>
  <ul class="dropdown-menu" role="menu"><li>
      <a href="../articles/cap-recap.html">Modeling variation in abundance using capture-recapture data</a>
    </li>
    <li>
      <a href="../articles/colext.html">Dynamic occupancy models in unmarked</a>
    </li>
    <li>
      <a href="../articles/contributing_to_unmarked.html">Contributing to unmarked: guide to adding a new model to `unmarked`</a>
    </li>
    <li>
      <a href="../articles/distsamp.html">Distance sampling analysis in unmarked</a>
    </li>
    <li>
      <a href="../articles/occuMulti.html">Multispecies occupancy models with occuMulti</a>
    </li>
    <li>
      <a href="../articles/powerAnalysis.html">Power Analysis in unmarked</a>
    </li>
    <li>
      <a href="../articles/simulate.html">Simulating datasets</a>
    </li>
    <li>
      <a href="../articles/spp-dist.html">Modeling and mapping species distributions</a>
    </li>
  </ul></li>
<li>
  <a href="../news/index.html">Changelog</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/hmecology/unmarked/" class="external-link">
    <span class="fab fa-github fa-lg"></span>

  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->



      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Class "unmarkedFitList"</h1>

    <div class="hidden name"><code>unmarkedFitList-class.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Class to hold multiple fitted models from one of
	<code>unmarked</code>'s fitting functions</p>
    </div>


    <div id="objects-from-the-class">
    <h2>Objects from the Class</h2>
    <p>Objects can be created by using the <code><a href="fitList.html">fitList</a></code> function.</p>
    </div>
    <div id="slots">
    <h2>Slots</h2>
    <p></p><dl><dt><code>fits</code>:</dt>
<dd><p>A <code>"list"</code> of models.</p></dd>


</dl></div>
    <div id="methods">
    <h2>Methods</h2>
    <p></p><dl><dt>coef</dt>
<dd><p><code>signature(object = "unmarkedFitList")</code>:
        Extract coefficients</p></dd>

    <dt>SE</dt>
<dd><p><code>signature(object = "unmarkedFitList")</code>:
        Extract standard errors</p></dd>

    <dt>modSel</dt>
<dd><p><code>signature(object = "unmarkedFitList")</code>:
		Model selection</p></dd>

    <dt>predict</dt>
<dd><p><code>signature(object = "unmarkedFitList")</code>:
		Model-averaged prediction</p></dd>

	
</dl></div>
    <div id="note">
    <h2>Note</h2>
    <p>Model-averaging regression coefficients is intentionally not implemented.</p>
    </div>
    <div id="see-also">
    <h2>See also</h2>
    <div class="dont-index"><p><code><a href="fitList.html">fitList</a></code>,
	<code><a href="unmarkedFit-class.html">unmarkedFit</a></code></p></div>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="fu">showClass</span><span class="op">(</span><span class="st">"unmarkedFitList"</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Class "unmarkedFitList" [package "unmarked"]</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Slots:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>            </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Name:  fits</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Class: list</span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/utils/data.html" class="external-link">data</a></span><span class="op">(</span><span class="va">linetran</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="op">(</span><span class="va">dbreaksLine</span> <span class="op">&lt;-</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">0</span>, <span class="fl">5</span>, <span class="fl">10</span>, <span class="fl">15</span>, <span class="fl">20</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> [1]  0  5 10 15 20</span>
<span class="r-in"><span><span class="va">lengths</span> <span class="op">&lt;-</span> <span class="va">linetran</span><span class="op">$</span><span class="va">Length</span> <span class="op">*</span> <span class="fl">1000</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">ltUMF</span> <span class="op">&lt;-</span> <span class="fu"><a href="https://rdrr.io/r/base/with.html" class="external-link">with</a></span><span class="op">(</span><span class="va">linetran</span>, <span class="op">{</span></span></span>
<span class="r-in"><span>  <span class="fu"><a href="unmarkedFrameDS.html">unmarkedFrameDS</a></span><span class="op">(</span>y <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/cbind.html" class="external-link">cbind</a></span><span class="op">(</span><span class="va">dc1</span>, <span class="va">dc2</span>, <span class="va">dc3</span>, <span class="va">dc4</span><span class="op">)</span>,</span></span>
<span class="r-in"><span>  siteCovs <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/data.frame.html" class="external-link">data.frame</a></span><span class="op">(</span><span class="va">Length</span>, <span class="va">area</span>, <span class="va">habitat</span><span class="op">)</span>, dist.breaks <span class="op">=</span> <span class="va">dbreaksLine</span>,</span></span>
<span class="r-in"><span>  tlength <span class="op">=</span> <span class="va">lengths</span>, survey <span class="op">=</span> <span class="st">"line"</span>, unitsIn <span class="op">=</span> <span class="st">"m"</span><span class="op">)</span></span></span>
<span class="r-in"><span>  <span class="op">}</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">fm1</span> <span class="op">&lt;-</span> <span class="fu"><a href="distsamp.html">distsamp</a></span><span class="op">(</span><span class="op">~</span> <span class="fl">1</span> <span class="op">~</span><span class="fl">1</span>, <span class="va">ltUMF</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">fm2</span> <span class="op">&lt;-</span> <span class="fu"><a href="distsamp.html">distsamp</a></span><span class="op">(</span><span class="op">~</span> <span class="va">area</span> <span class="op">~</span><span class="fl">1</span>, <span class="va">ltUMF</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">fm3</span> <span class="op">&lt;-</span> <span class="fu"><a href="distsamp.html">distsamp</a></span><span class="op">(</span> <span class="op">~</span> <span class="fl">1</span> <span class="op">~</span><span class="va">area</span>, <span class="va">ltUMF</span><span class="op">)</span></span></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">fl</span> <span class="op">&lt;-</span> <span class="fu"><a href="fitList.html">fitList</a></span><span class="op">(</span>Null<span class="op">=</span><span class="va">fm1</span>, A.<span class="op">=</span><span class="va">fm2</span>, .A<span class="op">=</span><span class="va">fm3</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">fl</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> An object of class "unmarkedFitList"</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Slot "fits":</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> $Null</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Call:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> distsamp(formula = ~1 ~ 1, data = ltUMF)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Density:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>  Estimate    SE     z P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>    -0.171 0.134 -1.28   0.201</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Detection:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>  Estimate    SE    z  P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>      2.39 0.127 18.7 2.46e-78</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> AIC: 164.7524 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> $A.</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Call:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> distsamp(formula = ~area ~ 1, data = ltUMF)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Density:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>  Estimate    SE     z P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>    -0.168 0.134 -1.25    0.21</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Detection:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>             Estimate     SE     z  P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> (Intercept)     3.00 0.5402  5.56 2.72e-08</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> area           -0.12 0.0955 -1.26 2.07e-01</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> AIC: 165.1845 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> $.A</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Call:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> distsamp(formula = ~1 ~ area, data = ltUMF)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Density:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>             Estimate     SE      z P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> (Intercept)   0.2364 0.5123  0.462   0.644</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> area         -0.0801 0.0979 -0.817   0.414</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Detection:</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>  Estimate    SE    z  P(&gt;|z|)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span>      2.39 0.127 18.7 2.47e-78</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> AIC: 166.0759 </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> </span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/stats/coef.html" class="external-link">coef</a></span><span class="op">(</span><span class="va">fl</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>        lam(Int)   p(Int)   p(area)   lam(area)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Null -0.1710554 2.386380        NA          NA</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> A.   -0.1678270 3.002507 -0.120364          NA</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> .A    0.2364320 2.386386        NA -0.08005895</span>
<span class="r-in"><span><span class="fu"><a href="SE-methods.html">SE</a></span><span class="op">(</span><span class="va">fl</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>       lam(Int)    p(Int)    p(area) lam(area)</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Null 0.1337819 0.1273598         NA        NA</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> A.   0.1340212 0.5401575 0.09548038        NA</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> .A   0.5122837 0.1273614         NA 0.0979427</span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">ms</span> <span class="op">&lt;-</span> <span class="fu"><a href="modSel.html">modSel</a></span><span class="op">(</span><span class="va">fl</span>, nullmod<span class="op">=</span><span class="st">"Null"</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">ms</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>      nPars    AIC delta AICwt cumltvWt   Rsq</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> Null     2 164.75  0.00  0.43     0.43 0.000</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> A.       3 165.18  0.43  0.35     0.78 0.122</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> .A       3 166.08  1.32  0.22     1.00 0.055</span>
<span class="r-in"><span></span></span>
<span class="r-in"><span></span></span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by <a href="https://chandlerlab.uga.edu/richard-chandler-phd/" class="external-link">Richard Chandler</a>, <a href="https://kenkellner.com" class="external-link">Ken Kellner</a>, Ian Fiske, David Miller, Andy Royle, Jeff Hostetler, Rebecca Hutchinson, Adam Smith, Lea Pautrel.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.1.0.</p>
</div>

      </footer></div>






  </body></html>

