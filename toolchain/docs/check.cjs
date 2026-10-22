// Copyright (C) The IMR developers. SPDX-License-Identifier: GPL-3.0-or-later
// Documentation-only checks. This file is never used by MATLAB simulations.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const http = require('node:http');
const {parse} = require('parse5');
const {chromium} = require('playwright');
const {mathjax} = require('mathjax-full/js/mathjax.js');
const {TeX} = require('mathjax-full/js/input/tex.js');
const {CHTML} = require('mathjax-full/js/output/chtml.js');
const {liteAdaptor} = require('mathjax-full/js/adaptors/liteAdaptor.js');
const {RegisterHTMLHandler} = require('mathjax-full/js/handlers/html.js');
const {AllPackages} = require('mathjax-full/js/input/tex/AllPackages.js');

const site = path.resolve(process.argv[2] || '../../site');
function files(directory) {
  return fs.readdirSync(directory, {withFileTypes: true}).flatMap(entry => {
    const file = path.join(directory, entry.name);
    return entry.isDirectory() ? files(file) : [file];
  });
}
function text(node) {
  return node.value || (node.childNodes || []).map(text).join('');
}
function walk(node, callback) {
  callback(node);
  for (const child of node.childNodes || []) walk(child, callback);
}

async function main() {
  const adaptor = liteAdaptor();
  RegisterHTMLHandler(adaptor);
  const math = mathjax.document('', {
    InputJax: new TeX({
      packages: AllPackages.filter(name => !['noerrors', 'noundefined'].includes(name)),
      formatError: (_, error) => { throw error; }
    }),
    OutputJax: new CHTML({fontURL: ''})
  });
  // Prove an invalid command fails rather than rendering a red fallback.
  assert.throws(() => math.convert('\\imrUndefinedCommand{z}'));
  let equations = 0;
  for (const file of files(site).filter(file => file.endsWith('.html'))) {
    walk(parse(fs.readFileSync(file, 'utf8')), node => {
      const attrs = Object.fromEntries((node.attrs || []).map(a => [a.name, a.value]));
      if (!(attrs.class || '').split(/\s+/).includes('arithmatex')) return;
      const tex = text(node).replace(/^\\[([]/, '').replace(/\\[)\]]$/, '');
      try {
        const rendered = adaptor.outerHTML(math.convert(tex, {display: node.tagName === 'div'}));
        assert(!rendered.includes('data-mjx-error'));
      } catch (error) {
        throw new Error(`${path.relative(site, file)}: ${tex}\n${error.message}`);
      }
      equations++;
    });
  }
  assert(equations > 0, 'No rendered equations were found.');
  console.log(`MathJax compiled ${equations} expressions; invalid-command guard passed.`);

  const mime = {'.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css',
    '.svg': 'image/svg+xml', '.woff': 'font/woff', '.woff2': 'font/woff2'};
  const server = http.createServer((request, response) => {
    const url = new URL(request.url, 'http://localhost');
    let file = path.resolve(site, '.' + decodeURIComponent(url.pathname));
    if (!(file === site || file.startsWith(site + path.sep))) {
      response.writeHead(403).end(); return;
    }
    if (fs.existsSync(file) && fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
    if (!fs.existsSync(file)) { response.writeHead(404).end(); return; }
    response.setHeader('Content-Type', mime[path.extname(file)] || 'application/octet-stream');
    fs.createReadStream(file).pipe(response);
  });
  let browser;
  try {
    await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
    browser = await chromium.launch({executablePath: process.env.IMR_CHROMIUM, args: ['--no-sandbox']});
    const page = await browser.newPage();
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    // Use the pinned MathJax 3 bundle for deterministic, offline rendering.
    const es5 = path.dirname(require.resolve('mathjax-full/es5/tex-mml-chtml.js'));
    await page.route('**/*', async route => {
      const url = new URL(route.request().url());
      if (url.hostname === '127.0.0.1') return route.continue();
      const marker = '/npm/mathjax@3/es5/';
      if (url.hostname === 'cdn.jsdelivr.net' && url.pathname.startsWith(marker)) {
        const file = path.resolve(es5, url.pathname.slice(marker.length));
        if (file.startsWith(es5 + path.sep) && fs.existsSync(file)) return route.fulfill({path: file});
      }
      return route.abort();
    });
    const base = `http://127.0.0.1:${server.address().port}`;
    const pages = ['theory/radial-equations/', 'model/heat-mass-transfer/', 'method/wall-boundaries/'];
    for (const [index, target] of pages.entries()) {
      if (index === 0) await page.goto(base + '/' + target);
      else {
        const link = page.locator(`a[href$="${target}"]`).first();
        await link.click();
        await page.waitForURL(base + '/' + target);
      }
      await page.waitForFunction(() => {
        const equations = [...document.querySelectorAll('.arithmatex')];
        return equations.length > 0 && equations.every(node => node.querySelector('mjx-container'));
      });
      assert.equal(await page.locator('mjx-merror, [data-mjx-error]').count(), 0);
      if (index === 0 && process.env.IMR_SCREENSHOT) await page.screenshot({path: process.env.IMR_SCREENSHOT});
    }
    assert.deepEqual(errors, []);
    console.log('Chromium rendered equations on initial load and two real page navigations.');
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
}
main().catch(error => { console.error(error); process.exitCode = 1; });
