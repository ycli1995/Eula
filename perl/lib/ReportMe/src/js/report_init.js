const config = window.REPORT_CONFIG;

function search(value, id) {
  $(id + " li").each(function () {
    var text = $(this).text();
    if (text.indexOf(value) == -1) {
      $(this).hide();
    } else {
      $(this).show();
    }
  });
}

function initNiceScroll(listSelector, containerSelector) {
  $(listSelector).niceScroll({
    cursoropacitymax: 0.5,
    cursorwidth: "8px",
  });
  $(containerSelector).niceScroll({
    cursoropacitymax: 0.5,
    cursorwidth: "8px",
  });
}

function onTabActivate(event) {
  var $tab = $(this);
  var $info = $("#nested-tabInfo2");
  var $name = $("span", $info);
  $name.text($tab.text());
  $info.show();
}

$(function () {
  $(".abbrTab").each(function () {
    var len = $(this).attr("data");
    var mess = $(this).text();
    $(this).attr("title", mess);
    if (mess.length > len - 1) {
      mess = mess.substring(0, len - 1);
      mess += '<span class="brief">...</span>';
    }
    $(this).html(mess);
  });

  $(".brief").click(function () {
    var info = $(this).parent().attr("title");
    messshow(info);
  });

  $(".func_table").DataTable({ ordering: false, scrollX: true });
});

$(document).ready(function () {
  $("#report_body").jumpto({
    innerWrapper: "section",
    firstLevel: "> h3",
    secondLevel: "> h5",
    offset: 0,
    anchorTopPadding: 90,
    animate: 600,
    showTitle: config.content_title,
    closeButton: false,
  });

  // init parent vertical tabs
  for (var i = 1; i < config.resp_vtabs_count; i++) {
    initNiceScroll("#resp-vtabs-list" + i, "#resp-vtabs-container" + i);
    $("#parentVerticalTab" + i).easyResponsiveTabs({
      type: "vertical",
      width: "auto",
      fit: true,
      closed: "accordion",
      tabidentify: "hor_" + i,
      activate: onTabActivate,
    });
  }

  for (var i = 1; i < config.resp_htabs_count; i++) {
    $("#parentHorizontalTab" + i).easyResponsiveTabs({
      type: "default",
      width: "auto",
      fit: true,
      closed: "accordion",
      tabidentify: "hor_" + i,
      activate: onTabActivate,
    });
  }

  // init child vertical tabs
  for (var i = 1; i < config.child_vtabs_count; i++) {
    initNiceScroll("#child-vtabs-list" + i, "#child-vtabs-container" + i);
    $("#ChildVerticalTab" + i).easyResponsiveTabs({
      type: "vertical",
      width: "auto",
      fit: true,
      tabidentify: "ver_" + i,
      activetab_bg: "#fff",
      inactive_bg: "#F5F5F5",
      active_border_color: "#c1c1c1",
      active_content_border_color: "#5AB1D0",
    });
  }
});

window.onload = function () {
  $(".jumpto-first").css("max-height", "calc(100vh - 150px)");
  $(".jumpto-first").niceScroll({
    cursoropacitymax: 0.5,
    cursorwidth: "6px",
    cursorborder: "0px",
  });
};

function messshow(info) {
  var str =
    '<div class="topBg"><div class="toppicBox"><div class="inBox" style="text-align:center;"><br><textarea class="pinfo" cols=56>' +
    info +
    '</textarea><br><br></div><img class="topClose" width="27" height="27" src="src/image/topclose.png" onclick="topclose()"></div></div>';
  $("body").append(str);
}

function topclose() {
  $(".topBg").remove();
}

function initTableHover(tableSelector) {
  $(tableSelector + " tr")
    .has("td")
    .each(function () {
      $(this).attr("onmouseover", 'this.style.backgroundColor = "#DDDDDD"');
      $(this).attr("onmouseout", 'this.style.backgroundColor = "#FFFFFF"');
    });
}

initTableHover(".hl_table");
initTableHover(".func_table");

$(".logo").click(function () {
  $(window).scrollTop(0);
});

var width = $(".func_table:eq(0)").parents("section").width();
if (width < 200) {
  width = 1000;
}
$(".func_table caption").css("width", width);
$(".func_table")
  .parents(".dataTables_scrollBody")
  .scroll(function () {
    $(this).parent().find("caption").css("margin-left", $(this).scrollLeft());
  });
